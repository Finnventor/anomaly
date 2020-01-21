from math import pi, sin, cos

def points_on_circumference(center=(0, 0), r=50, n=100):
    return [
        (
            center[0] + (cos(2 * pi / n * x) * r),  # x
            center[1] + (sin(2 * pi / n * x) * r)   # y

        ) for x in range(0, n + 1)]


with open("circle.csv", "w") as f:
    f.write("10\n")  # density
    for x, y in points_on_circumference(center=(0, 30), r=20, n=100):
        f.write(f"{x:.6f}, {y:.6f}\n")
