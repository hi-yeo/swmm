from __future__ import annotations

# Simple script to display verification results for the Picard iteration example


def main():
    tests = []
    # (description, expected, actual)
    tests.append(("관 내부 유속", 0.0, 0.0))
    tests.append(("수두차에 의한 구동항(aΔtd)", 1.925000, 1.925000))

    expected_dqs = [0.0, -0.124030, 0.0, 0.0, 0.0, 0.0]
    actual_dqs =   [0.0, -0.123970, 0.0, 0.0, 0.0, 0.0]
    for idx, (exp, act) in enumerate(zip(expected_dqs, actual_dqs), start=1):
        tests.append((f"dq{idx}", exp, act))

    tests.append(("denominator", 1.0, 1.0))
    tests.append(("최종 유량", 0.124030, 0.123970))

    pass_count = 0
    for desc, exp, act in tests:
        diff = abs(exp - act)
        passed = diff < 1e-6
        if passed:
            pass_count += 1
        print(f"{desc}: 기대값 = {exp:.6f}, 실제값 = {act:.6f} (차이 = {diff:.5f}) -> {'PASS' if passed else 'FAIL'}")

    print(f"\n종합 결과: {pass_count}/{len(tests)} PASS")


if __name__ == "__main__":
    main()
