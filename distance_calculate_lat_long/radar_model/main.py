from input_handler import read_point
from radar_calculations import calculate_radar_parameters
from intermediate_points import calculate_intermediate_points


def main():
    # Step 1: Input
    P1 = read_point("P1")
    P2 = read_point("P2")

    n = int(input("\nNumber of intermediate points: "))

    # Step 2: Radar calculations
    result = calculate_radar_parameters(P1, P2)

    print("\n========== RADAR RESULTS ==========")
    print(f"Slant Range: {result['slant_range']:.3f} m")
    print(f"Azimuth: {result['azimuth']:.3f} degrees")
    print(f"Elevation: {result['elevation']:.3f} degrees")

    # Step 3: Intermediate coordinates
    points = calculate_intermediate_points(P1, P2, n)

    print("\n======= INTERMEDIATE POINTS =======")

    for i, point in enumerate(points, start=1):
        print(
            f"Point {i}: "
            f"Lat={point['latitude']:.6f}, "
            f"Lon={point['longitude']:.6f}, "
            f"Height={point['height']:.3f} m"
        )


if __name__ == "__main__":
    main()