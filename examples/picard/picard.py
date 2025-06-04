import ctypes
import math
import os
import subprocess

GRAVITY = 32.2
PHI = 1.486


def area(width, y):
    return width * y

def hyd_radius(width, y):
    return width * y / (width + 2 * y)

def compute_step(q_last, y1, y2, h1, h2, q_old, a_old, dt, length, width, n, omega):
    a1 = area(width, y1)
    a2 = area(width, y2)
    a_mid = 0.5 * (a1 + a2)
    r1 = hyd_radius(width, y1)
    r_mid = hyd_radius(width, 0.5 * (y1 + y2))

    rho = 1.0
    a_wtd = a1 + (a_mid - a1) * rho
    r_wtd = r1 + (r_mid - r1) * rho

    v = q_last / a_mid
    rough_factor = GRAVITY * (n / PHI) ** 2

    dq1 = dt * rough_factor / (r_wtd ** 1.33333) * abs(v)
    dq2 = dt * GRAVITY * a_wtd * (h2 - h1) / length
    denom = 1.0 + dq1
    q = (q_old - dq2) / denom

    q = (1.0 - omega) * q_last + omega * q
    return q

def picard_flow():
    width = 3.0
    length = 200.0
    n = 0.013
    dt = 1.0
    omega = 0.5
    y1 = 5.0
    y2 = 4.9
    h1 = 5.0
    h2 = 4.9
    q_old = 0.0
    a_old = area(width, (y1 + y2) / 2.0)

    q_last = q_old
    for _ in range(20):
        q_new = compute_step(q_last, y1, y2, h1, h2, q_old, a_old, dt, length, width, n, omega)
        if abs(q_new - q_last) < 1e-6:
            break
        q_last = q_new
    return q_last


def run_c_binary():
    c_path = os.path.join(os.path.dirname(__file__), "picard.c")
    exe_path = os.path.join(os.path.dirname(__file__), "picard_c")
    subprocess.check_call(["gcc", c_path, "-lm", "-o", exe_path])
    output = subprocess.check_output([exe_path])
    return float(output.decode().strip())

if __name__ == "__main__":
    py_q = picard_flow()
    c_q = run_c_binary()
    print(f"Python result: {py_q:.6f}")
    print(f"C result: {c_q:.6f}")
    print(f"Difference: {abs(py_q - c_q):.6e}")
