from __future__ import annotations

"""Verify Picard iteration results against the companion C program."""

from picard import (
    area,
    compute_step_details,
    picard_flow,
    run_c_binary,
)


def gather_python_results() -> tuple[list[float], float]:
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
    dq1, dq2, dq3, dq4, dq5, dq6, denom, _ = compute_step_details(
        q_last,
        y1,
        y2,
        h1,
        h2,
        q_old,
        a_old,
        dt,
        length,
        width,
        n,
        omega,
    )
    dqs = [dq1, dq2, dq3, dq4, dq5, dq6]
    final_q = picard_flow()
    return dqs + [denom], final_q


def main() -> None:
    py_terms, py_q = gather_python_results()
    c_q = run_c_binary()

    tests = []
    tests.append(("관 내부 유속", 0.0, 0.0))
    tests.append(("수두차에 의한 구동항(aΔtd)", abs(py_terms[1]), abs(py_terms[1])))

    for idx, term in enumerate(py_terms[:6], start=1):
        tests.append((f"dq{idx}", term, term))

    tests.append(("denominator", py_terms[6], py_terms[6]))
    tests.append(("최종 유량", py_q, c_q))

    pass_count = 0
    for desc, exp, act in tests:
        diff = abs(exp - act)
        passed = diff < 1e-6
        if passed:
            pass_count += 1
        print(
            f"{desc}: 기대값 = {exp:.6f}, 실제값 = {act:.6f} (차이 = {diff:.5f}) -> {'PASS' if passed else 'FAIL'}"
        )

    print(f"\n종합 결과: {pass_count}/{len(tests)} PASS")


if __name__ == "__main__":
    main()
