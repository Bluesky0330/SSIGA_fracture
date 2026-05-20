import numpy as np
from scipy import special

def calculate_elliptic_integrals_2nd(phi_deg, k):
    """
    第二種楕円積分を計算する関数
    
    Parameters:
        phi_deg (float): 振幅（度数法）
        k (float): 母数 (modulus, 0 <= k < 1)
        
    Returns:
        tuple: (不完全楕円積分の値, 完全楕円積分の値)
    """
    # 角度をラジアンに変換
    phi_rad = np.radians(phi_deg)
    
    # SciPyでは m = k^2 を引数として渡す必要がある
    m = k**2
    
    # 第二種不完全楕円積分 E(phi, m)
    # scipy.special.ellipeinc(phi, m)
    inc_val = special.ellipeinc(phi_rad, m)
    
    # 第二種完全楕円積分 E(m)
    # scipy.special.ellipe(m)
    com_val = special.ellipe(m)
    
    return inc_val, com_val

# --- 使用例 ---

# 設定値
k_val = 0.979795897113271   # 母数
phi_val = 45  # 角度（度）

# 計算実行
inc_result, com_result = calculate_elliptic_integrals_2nd(phi_val, k_val)

print(f"母数 k = {k_val}, 角度 phi = {phi_val}° の場合:")
print("-" * 40)
print(f"第二種 不完全楕円積分 E(phi, k): {inc_result:.20f}")
print(f"第二種 完全楕円積分 E(k):       {com_result:.20f}")

# 参考: 楕円の周長計算への応用 (長半径 a, 離心率 e=k)
# 周長 L = 4 * a * E(e)
a = 10.0
L = 4 * a * com_result
print("-" * 40)
print(f"長半径 a={a}, 離心率 e={k_val} の楕円の周長: {L:.6f}")