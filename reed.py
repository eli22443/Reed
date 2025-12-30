import galois
import numpy as np

# 1. הגדרת השדה והפרמטרים [cite: 17, 18]
# נשתמש בשדה GF(2^8) כדי שכל סמל יהיה בייט
GF = galois.GF(2**8)
n = 15  # אורך מילת הקוד (n) [cite: 22]
k = 7   # אורך הודעת המקור (k) [cite: 22]
alpha = GF(list(range(n)))  # נקודות הערכה [cite: 18]

def encode_rs(message_symbols):
    """קידוד הודעה לפולינום והערכתו בנקודות alpha [cite: 20, 21]"""
    # יצירת פולינום M(x) מהודעת המקור 
    poly = galois.Poly(message_symbols[::-1], field=GF)
    # הערכה בנקודות ליצירת Codeword [cite: 21]
    return poly(alpha)

def sudan_list_decoding(received_y, n, k):
    dy = 1
    dx = k + 2 
    
    # בניית המטריצה (אינטרפולציה) [cite: 30]
    A = []
    for i in range(n):
        row = []
        for j in range(dy + 1):
            for r in range(dx + 1):
                row.append((alpha[i]**r) * (received_y[i]**j))
        A.append(row)
    
    A = GF(A)
    null_space = A.null_space()
    # ... (המשך אחרי null_space = A.null_space())

    if len(null_space) == 0:
        return []

    # לוקחים את פתרון המקדמים הראשון שמצאנו [cite: 31]
    coeffs = null_space[0]
    
    # חילוץ המקדמים של A(x) ו-B(x) לפי סדר המילוי במטריצה שלכם
    # במטריצה שלכם (שורות 26-27): r רץ בתוך j
    # לכן dx+1 המקדמים הראשונים שייכים ל-j=0 (כלומר A(x))
    # והבאים שייכים ל-j=1 (כלומר B(x))
    num_x = dx + 1
    coeffs_A = coeffs[0 : num_x]
    coeffs_B = coeffs[num_x : 2 * num_x]

    # יצירת הפולינומים [cite: 20]
    A_x = galois.Poly(coeffs_A[::-1], field=GF)
    B_x = galois.Poly(coeffs_B[::-1], field=GF)

    # מציאת הודעת המקור p(X) = -A(X) / B(X) 
    try:
        # חילוק פולינומים
        p_x, remainder = divmod(-A_x, B_x)
        
        # אם החילוק מדויק והדרגה מתאימה, מצאנו מועמד! [cite: 33]
        if remainder == 0 and p_x.degree <= k-1:
            decoded_msg = p_x.coeffs(size=k)[::-1]
            return [decoded_msg]
    except:
        pass

    return []

# --- הרצה ובדיקה ---

# א. יצירת הודעה מקורית [cite: 19]
original_msg = GF([10, 20, 30, 40, 50, 60, 70])
print(f"Original Message: {original_msg}")

# ב. קידום [cite: 21]
codeword = encode_rs(original_msg)

# ג. הוספת רעש (מעבר ליכולת של פענוח ייחודי) [cite: 11, 14]
# ייחודי יכול לתקן עד (15-7)/2 = 4 שגיאות [cite: 8]
# אנחנו נזריק 6 שגיאות כדי להראות את כוחו של List Decoding [cite: 14]
received = codeword.copy()
error_indices = [0, 1, 2, 3, 4, 5]
for idx in error_indices:
    received[idx] = GF.Random()

print(f"Received with {len(error_indices)} errors.")

# ד. פענוח רשימה [cite: 9]
candidates = sudan_list_decoding(received, n, k)
print(f"Success! Found {len(candidates)} potential messages in the list.")
print(candidates)