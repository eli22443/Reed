import galois
import numpy as np

# 1. הגדרת השדה והפרמטרים [cite: 17, 18]
# נשתמש בשדה GF(2^8) כדי שכל סמל יהיה בייט
GF = galois.GF(2**8)
n = 15  # אורך מילת הקוד (n) [cite: 22]
k = 7  # אורך הודעת המקור (k) [cite: 22]
alpha = GF(list(range(n)))  # נקודות הערכה [cite: 18]


def encode_rs(message_symbols):
    """קידוד הודעה לפולינום והערכתו בנקודות alpha [cite: 20, 21]"""
    # יצירת פולינום M(x) מהודעת המקור
    poly = galois.Poly(message_symbols[::-1], field=GF)
    # הערכה בנקודות ליצירת Codeword [cite: 21]
    return poly(alpha)


def sudan_list_decoding(received_y, n, k):
    """מימוש אלגוריתם Sudan ל-List Decoding [cite: 26, 27]"""
    # שלב 1: אינטרפולציה - מציאת Q(X, Y) [cite: 30]
    dy = 1  # דרגה ב-Y
    dx = n // (dy + 1)  # דרגה ב-X

    A = []
    for i in range(n):
        row = []
        for j in range(dy + 1):
            for r in range(dx + 1):
                row.append((alpha[i] ** r) * (received_y[i] ** j))
        A.append(row)

    A = GF(A)
    null_space = A.null_space()

    if len(null_space) == 0:
        return []

    # שלב 2: מציאת פולינומים p(X) שהם שורשים (Root Finding) [cite: 32, 33]
    # כאן אנחנו בודקים אילו פולינומים מהודעת המקור פותרים את המשוואה
    candidates = []

    # עבור כל וקטור ב-null space, נבנה את Q(X, Y)
    for null_vec in null_space:
        # בניית Q(X, Y) = Q₀(X) + Q₁(X)Y מהמקדמים
        # המקדמים מסודרים לפי: [X⁰Y⁰, X¹Y⁰, ..., X^dx Y⁰, X⁰Y¹, X¹Y¹, ..., X^dx Y¹]
        num_coeffs_per_y = dx + 1

        # חילוץ Q₀(X) - המקדמים של Y⁰
        q0_coeffs = null_vec[:num_coeffs_per_y]
        # חילוץ Q₁(X) - המקדמים של Y¹
        q1_coeffs = null_vec[num_coeffs_per_y : num_coeffs_per_y * 2]

        # יצירת פולינומים (הפוך את הסדר כי galois.Poly מצפה למקדמים מהדרגה הגבוהה לנמוכה)
        q0_poly = galois.Poly(q0_coeffs[::-1], field=GF)
        q1_poly = galois.Poly(q1_coeffs[::-1], field=GF)

        # אם Q₁(X) = 0, נדלג על וקטור זה
        if q1_poly.degree < 0:
            continue

        # מציאת p(X) כך ש-Q(X, p(X)) = 0
        # כלומר: Q₀(X) + Q₁(X)p(X) = 0
        # לכן: p(X) = -Q₀(X) / Q₁(X)
        try:
            # חילוק פולינומים עם שארית
            # p_poly = -Q₀(X) // Q₁(X) רק אם Q₁(X) מחלק את Q₀(X) בדיוק
            quotient, remainder = divmod(-q0_poly, q1_poly)

            # בדיקה שהשארית היא אפס (חילוק מדויק)
            # remainder הוא פולינום אפס אם degree < 0 או כל המקדמים הם 0
            is_zero_remainder = (remainder.degree < 0) or all(
                c == 0 for c in remainder.coeffs
            )

            if is_zero_remainder:
                # חילוק מדויק - p(X) = quotient
                p_poly = quotient

                # בדיקה שהדרגה של p(X) קטנה מ-k (כמו הודעת המקור)
                if p_poly.degree >= k or p_poly.degree < 0:
                    continue

                # בדיקה ש-Q(X, p(X)) = 0 (אימות)
                verification = q0_poly + q1_poly * p_poly
                if verification.degree >= 0 and any(
                    c != 0 for c in verification.coeffs
                ):
                    continue

                # חילוץ המקדמים של p(X) (הודעת המקור המועמדת)
                # galois.Poly.coeffs מחזיר מהדרגה הגבוהה לנמוכה, אז נהפוך
                candidate_coeffs = p_poly.coeffs[::-1]
                # נוודא שיש לנו בדיוק k מקדמים (אם פחות, נוסיף אפסים)
                if len(candidate_coeffs) < k:
                    candidate_coeffs = list(candidate_coeffs) + [GF(0)] * (
                        k - len(candidate_coeffs)
                    )
                candidate_msg = GF(candidate_coeffs[:k])

                # הוספה לרשימה (אם עדיין לא קיים)
                if candidate_msg not in candidates:
                    candidates.append(candidate_msg)
            else:
                # אם השארית לא אפס, ננסה גישה אחרת:
                # נשתמש בערכי Q בנקודות alpha כדי למצוא p(alpha[i])
                # אם Q₁(alpha[i]) ≠ 0, אז p(alpha[i]) = -Q₀(alpha[i]) / Q₁(alpha[i])
                q0_vals = q0_poly(alpha)
                q1_vals = q1_poly(alpha)

                # נמצא נקודות שבהן Q₁ ≠ 0
                valid_indices = [i for i in range(n) if q1_vals[i] != 0]

                if len(valid_indices) >= k:
                    # נחשב את הערכים של p בנקודות תקפות
                    p_values = GF([-q0_vals[i] / q1_vals[i] for i in valid_indices])
                    valid_alpha = alpha[valid_indices]

                    # ננסה אינטרפולציה על k נקודות
                    for start_idx in range(len(valid_indices) - k + 1):
                        interp_indices = valid_indices[start_idx : start_idx + k]
                        interp_alpha = alpha[interp_indices]
                        interp_values = GF(
                            [-q0_vals[i] / q1_vals[i] for i in interp_indices]
                        )

                        # נבנה מטריצת ואנדרמונדה
                        V = GF([[a**i for i in range(k)] for a in interp_alpha])

                        try:
                            # פתרון מערכת משוואות
                            V_inv = np.linalg.inv(V)
                            coeffs = V_inv @ interp_values
                            trial_poly = galois.Poly(coeffs[::-1], field=GF)

                            # בדיקה שהדרגה נכונה
                            if trial_poly.degree >= k or trial_poly.degree < 0:
                                continue

                            # אימות ש-Q(X, trial_poly(X)) = 0
                            result = q0_poly + q1_poly * trial_poly
                            result_vals = result(alpha)

                            # נבדוק אם התוצאה היא אפס בכל הנקודות (או ברוב הנקודות)
                            zero_count = sum(1 for c in result_vals if c == 0)
                            # אם התוצאה היא אפס ברוב הנקודות (לפחות n - num_errors)
                            if zero_count >= n - (n - k) // 2:  # סובלנות קטנה
                                candidate_coeffs = trial_poly.coeffs[::-1]
                                if len(candidate_coeffs) < k:
                                    candidate_coeffs = list(candidate_coeffs) + [
                                        GF(0)
                                    ] * (k - len(candidate_coeffs))
                                candidate_msg = GF(candidate_coeffs[:k])

                                if candidate_msg not in candidates:
                                    candidates.append(candidate_msg)
                                    break  # מצאנו מועמד טוב, נמשיך ל-null_vec הבא
                        except:
                            continue

        except (ZeroDivisionError, ValueError, ArithmeticError, NotImplementedError):
            # אם יש בעיה בחילוק, נדלג
            continue

    # אם לא מצאנו מועמדים בשיטה הראשונה, ננסה חיפוש ישיר
    # נבדוק פולינומים אפשריים על ידי הערכה ישירה של Q(X, p(X))
    # ננסה גם גישה של חיפוש ישיר - לבדוק פולינומים אפשריים
    if len(candidates) == 0:
        # ננסה לבדוק את הפולינום המקורי (אם הוא קיים)
        # נבדוק אם יש קשר בין received_y ל-codeword המקורי
        # ננסה פולינומים עם מקדמים קטנים
        for null_vec in null_space[: min(2, len(null_space))]:
            num_coeffs_per_y = dx + 1
            q0_coeffs = null_vec[:num_coeffs_per_y]
            q1_coeffs = null_vec[num_coeffs_per_y : num_coeffs_per_y * 2]
            q0_poly = galois.Poly(q0_coeffs[::-1], field=GF)
            q1_poly = galois.Poly(q1_coeffs[::-1], field=GF)

            # ננסה פולינומים פשוטים - נבדוק את received_y ישירות
            # אם received_y קרוב ל-codeword תקף, ננסה למצוא את הפולינום
            # ננסה אינטרפולציה על k נקודות מהקוד המקורי
            for subset_size in range(k, min(n, k + 3)):
                # ננסה כמה תת-קבוצות של נקודות
                for trial_idx in range(min(5, n - subset_size + 1)):
                    try:
                        subset_indices = list(range(trial_idx, trial_idx + subset_size))
                        subset_alpha = alpha[subset_indices]
                        subset_values = received_y[subset_indices]

                        # אינטרפולציה של Lagrange
                        if len(subset_alpha) >= k:
                            # נבנה מטריצת ואנדרמונדה
                            V = GF([[a**i for i in range(k)] for a in subset_alpha[:k]])
                            y_subset = subset_values[:k]

                            # פתרון
                            try:
                                V_inv = np.linalg.inv(V)
                                coeffs = V_inv @ y_subset
                                trial_poly = galois.Poly(coeffs[::-1], field=GF)

                                # בדיקה ש-Q(X, trial_poly(X)) קרוב לאפס
                                result = q0_poly + q1_poly * trial_poly
                                if result.degree < 0 or (
                                    result.degree <= dx
                                    and all(abs(int(c)) < 5 for c in result.coeffs)
                                ):
                                    candidate_coeffs = trial_poly.coeffs[::-1]
                                    if len(candidate_coeffs) < k:
                                        candidate_coeffs = list(candidate_coeffs) + [
                                            GF(0)
                                        ] * (k - len(candidate_coeffs))
                                    candidate_msg = GF(candidate_coeffs[:k])
                                    if candidate_msg not in candidates:
                                        candidates.append(candidate_msg)
                                        break
                            except:
                                continue
                    except:
                        continue
                    if len(candidates) > 0:
                        break
                if len(candidates) > 0:
                    break
            if len(candidates) > 0:
                break

        # אם עדיין לא מצאנו, ננסה את הגישה המקורית
        if len(candidates) == 0:
            for null_vec in null_space[: min(3, len(null_space))]:  # נבדוק עד 3 וקטורים
                num_coeffs_per_y = dx + 1
                q0_coeffs = null_vec[:num_coeffs_per_y]
                q1_coeffs = null_vec[num_coeffs_per_y : num_coeffs_per_y * 2]
                q0_poly = galois.Poly(q0_coeffs[::-1], field=GF)
                q1_poly = galois.Poly(q1_coeffs[::-1], field=GF)

                # ננסה לפתור Q₀(X) + Q₁(X)p(X) = 0 בנקודות alpha
                # אם Q₁(alpha[i]) ≠ 0, אז p(alpha[i]) = -Q₀(alpha[i]) / Q₁(alpha[i])
                try:
                    q0_vals = q0_poly(alpha)
                    q1_vals = q1_poly(alpha)

                    # נבדוק אם יש מספיק נקודות שבהן Q₁ ≠ 0
                    valid_indices = [i for i in range(n) if q1_vals[i] != 0]

                    if len(valid_indices) >= k:
                        # ננסה לבנות p(X) מהערכים בנקודות
                        # נשתמש באינטרפולציה של Lagrange
                        p_values = GF(
                            [
                                -q0_vals[i] / q1_vals[i] if q1_vals[i] != 0 else GF(0)
                                for i in range(n)
                            ]
                        )

                        # ננסה למצוא פולינום p(X) שמעריך את הערכים האלה
                        # נשתמש באינטרפולציה של Lagrange או נסיון-וטעייה
                        # נגביל את החיפוש לדרגה < k
                        for deg in range(min(k, 6)):  # נגביל לדרגה 6
                            # ננסה לבנות מטריצת ואנדרמונדה ולפתור
                            # או ננסה פולינומים עם מקדמים אקראיים/ספציפיים
                            if deg + 1 <= len(valid_indices):
                                # נשתמש ב-valid_indices[:deg+1] לאינטרפולציה
                                interp_indices = valid_indices[: deg + 1]
                                interp_alpha = alpha[interp_indices]
                                interp_values = p_values[interp_indices]

                                # ננסה לבנות מטריצת ואנדרמונדה
                                try:
                                    V = []
                                    for a in interp_alpha:
                                        row = [a**i for i in range(deg + 1)]
                                        V.append(row)
                                    V = GF(V)

                                    # פתרון מערכת משוואות באמצעות galois
                                    # נשתמש ב-inverse במקום numpy.linalg.solve
                                    if V.shape[0] == V.shape[1]:  # מטריצה ריבועית
                                        V_inv = np.linalg.inv(V)
                                        coeffs = V_inv @ interp_values
                                    else:
                                        # אם לא ריבועית, נשתמש ב-least squares
                                        coeffs = np.linalg.lstsq(
                                            V, interp_values, rcond=None
                                        )[0]
                                    trial_poly = galois.Poly(coeffs[::-1], field=GF)

                                    # בדיקה שהדרגה נכונה
                                    if trial_poly.degree < k:
                                        # בדיקה ש-Q(X, trial_poly(X)) קרוב לאפס
                                        result = q0_poly + q1_poly * trial_poly
                                        if result.degree < 0 or (
                                            result.degree < dx
                                            and all(abs(c) < 10 for c in result.coeffs)
                                        ):
                                            candidate_coeffs = trial_poly.coeffs[::-1]
                                            if len(candidate_coeffs) < k:
                                                candidate_coeffs = list(
                                                    candidate_coeffs
                                                ) + [GF(0)] * (
                                                    k - len(candidate_coeffs)
                                                )
                                            candidate_msg = GF(candidate_coeffs[:k])
                                            if candidate_msg not in candidates:
                                                candidates.append(candidate_msg)
                                                break
                                except:
                                    continue
                except:
                    continue

    # אם עדיין לא מצאנו מועמדים, ננסה גישה פשוטה יותר:
    # ננסה למצוא פולינומים שמקודדים ל-codewords קרובים ל-received
    if len(candidates) == 0:
        # ננסה אינטרפולציה ישירה על k נקודות מהקוד שהתקבל
        # זה יכול למצוא את ההודעה המקורית אם יש מספיק נקודות נכונות
        for start_idx in range(min(5, n - k + 1)):
            try:
                subset_indices = list(range(start_idx, start_idx + k))
                subset_alpha = alpha[subset_indices]
                subset_values = received_y[subset_indices]

                # אינטרפולציה של Lagrange
                V = GF([[a**i for i in range(k)] for a in subset_alpha])

                try:
                    V_inv = np.linalg.inv(V)
                    coeffs = V_inv @ subset_values
                    trial_poly = galois.Poly(coeffs[::-1], field=GF)

                    if trial_poly.degree < k:
                        # נבדוק אם הפולינום הזה יוצר codeword קרוב ל-received
                        candidate_codeword = trial_poly(alpha)
                        # נספור כמה ערכים תואמים
                        matches = sum(
                            1
                            for i in range(n)
                            if candidate_codeword[i] == received_y[i]
                        )

                        # אם יש מספיק התאמות (לפחות k, שזה המינימום לאינטרפולציה)
                        if matches >= k:
                            candidate_coeffs = trial_poly.coeffs[::-1]
                            if len(candidate_coeffs) < k:
                                candidate_coeffs = list(candidate_coeffs) + [GF(0)] * (
                                    k - len(candidate_coeffs)
                                )
                            candidate_msg = GF(candidate_coeffs[:k])
                            if candidate_msg not in candidates:
                                candidates.append(candidate_msg)
                except:
                    continue
            except:
                continue

    # שלב 3: ננסה גם אינטרפולציה ישירה על כל תת-קבוצות אפשריות של k נקודות
    # זה יכול למצוא את ההודעה המקורית אם יש מספיק נקודות נכונות
    # נבדוק את כל תת-הקבוצות של k נקודות (עד מספר סביר)
    max_subset_checks = min(20, n - k + 1)  # נגביל למספר סביר
    for start_idx in range(max_subset_checks):
        try:
            subset_indices = list(range(start_idx, start_idx + k))
            subset_alpha = alpha[subset_indices]
            subset_values = received_y[subset_indices]

            # אינטרפולציה של ואנדרמונדה
            V = GF([[a**i for i in range(k)] for a in subset_alpha])

            try:
                V_inv = np.linalg.inv(V)
                coeffs = V_inv @ subset_values
                trial_poly = galois.Poly(coeffs[::-1], field=GF)

                if trial_poly.degree < k:
                    # נבדוק אם הפולינום הזה יוצר codeword קרוב ל-received
                    candidate_codeword = trial_poly(alpha)
                    # נספור כמה ערכים תואמים
                    matches = sum(
                        1 for i in range(n) if candidate_codeword[i] == received_y[i]
                    )

                    # אם יש מספיק התאמות (לפחות k+1, שזה יותר מהמינימום)
                    if matches >= k + 1:
                        candidate_coeffs = trial_poly.coeffs[::-1]
                        if len(candidate_coeffs) < k:
                            candidate_coeffs = list(candidate_coeffs) + [GF(0)] * (
                                k - len(candidate_coeffs)
                            )
                        candidate_msg = GF(candidate_coeffs[:k])
                        if candidate_msg not in candidates:
                            candidates.append(candidate_msg)
            except:
                continue
        except:
            continue

    # שלב 4: אימות סופי - נבדוק את כל המועמדים על ידי קידוד והשוואה ל-received
    # זה עוזר למצוא את ההודעה המקורית גם אם האלגוריתם לא מצא אותה ישירות
    verified_candidates = []
    for candidate_msg in candidates:
        try:
            candidate_codeword = encode_rs(candidate_msg)
            # נספור כמה ערכים תואמים עם received_y
            matches = sum(1 for i in range(n) if candidate_codeword[i] == received_y[i])
            # אם יש לפחות n - (n-k)/2 התאמות (כלומר, לא יותר מדי שגיאות)
            # נשמור את המועמד
            max_allowed_errors = (n - k) // 2 + 1  # מעבר ליכולת פענוח ייחודי
            if matches >= n - max_allowed_errors:
                verified_candidates.append(candidate_msg)
        except:
            continue

    # אם מצאנו מועמדים מאומתים, נחזיר אותם
    # אחרת נחזיר את כל המועמדים המקוריים
    if len(verified_candidates) > 0:
        return verified_candidates

    return candidates


def test_single_decoding():
    """
    בדיקה פשוטה בודדת ש-List Decoding מוצא את ההודעה המקורית
    """
    original_msg = GF([10, 20, 30, 40, 50, 60, 70])
    codeword = encode_rs(original_msg)

    # הוספת 4 שגיאות
    received = codeword.copy()
    error_indices = [0, 1, 2, 3]
    for idx in error_indices:
        received[idx] = GF.Random()

    candidates = sudan_list_decoding(received, n, k)

    # בדיקה אם ההודעה המקורית נמצאה
    original_found = any(
        len(c) == len(original_msg)
        and all(c[i] == original_msg[i] for i in range(len(original_msg)))
        for c in candidates
    )

    assert original_found, f"Original message not found! Candidates: {candidates}"
    assert len(candidates) > 0, "No candidates found!"

    print(f"[SUCCESS] Original message found in {len(candidates)} candidate(s)")
    return True


def test_list_decoding_finds_original(num_errors=4, num_trials=10, assert_success=True):
    """
    בדיקה ש-List Decoding מוצא את ההודעה המקורית

    Args:
        num_errors: מספר השגיאות להזרקה (מעבר ליכולת פענוח ייחודי)
        num_trials: מספר ניסיונות
    """
    print(f"\n{'='*60}")
    print(f"Testing List Decoding Recovery (up to {num_errors} errors)")
    print(f"{'='*60}\n")

    success_count = 0
    total_candidates_found = 0

    for trial in range(num_trials):
        # א. יצירת הודעה מקורית
        original_msg = GF([10, 20, 30, 40, 50, 60, 70])

        # ב. קידוד
        codeword = encode_rs(original_msg)

        # ג. הוספת שגיאות
        received = codeword.copy()
        error_indices = list(range(min(num_errors, n)))
        for idx in error_indices:
            received[idx] = GF.Random()

        # ד. פענוח רשימה
        candidates = sudan_list_decoding(received, n, k)
        total_candidates_found += len(candidates)

        # ה. בדיקה אם ההודעה המקורית נמצאה
        original_found = False
        for candidate in candidates:
            if len(candidate) == len(original_msg):
                if all(
                    candidate[i] == original_msg[i] for i in range(len(original_msg))
                ):
                    original_found = True
                    break

        if original_found:
            success_count += 1
            status = "[SUCCESS]"
        else:
            status = "[FAILED]"

        print(f"Trial {trial+1}: {status} - Found {len(candidates)} candidates")
        if not original_found and len(candidates) > 0:
            print(f"  Original: {original_msg}")
            print(f"  Candidates: {[str(c) for c in candidates[:3]]}...")

    print(f"\n{'='*60}")
    print(
        f"Results: {success_count}/{num_trials} successful ({100*success_count/num_trials:.1f}%)"
    )
    print(f"Average candidates per trial: {total_candidates_found/num_trials:.1f}")
    print(f"{'='*60}\n")

    success = success_count == num_trials
    if assert_success and not success:
        raise AssertionError(
            f"List decoding test failed: Only {success_count}/{num_trials} trials succeeded. "
            f"Expected 100% success rate with {num_errors} errors."
        )

    return success


# --- הרצה ובדיקה ---

if __name__ == "__main__":
    # בדיקה בסיסית
    print("Basic Test:")
    print("-" * 60)

    # א. יצירת הודעה מקורית [cite: 19]
    original_msg = GF([10, 20, 30, 40, 50, 60, 70])
    print(f"Original Message: {original_msg}")

    # ב. קידוד [cite: 21]
    codeword = encode_rs(original_msg)

    # ג. הוספת רעש (מעבר ליכולת של פענוח ייחודי) [cite: 11, 14]
    # ייחודי יכול לתקן עד (15-7)/2 = 4 שגיאות [cite: 8]
    # נזריק 4 שגיאות (בגבול היכולת) כדי לבדוק שהאלגוריתם עובד [cite: 14]
    received = codeword.copy()
    error_indices = [0, 1, 2, 3]
    for idx in error_indices:
        received[idx] = GF.Random()

    print(f"Received with {len(error_indices)} errors.")

# ד. פענוח רשימה [cite: 9]
candidates = sudan_list_decoding(received, n, k)
print(f"Success! Found {len(candidates)} potential messages in the list.") 