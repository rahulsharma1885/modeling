from geographiclib.geodesic import Geodesic


def calculate_intermediate_points(P1, P2, n):
    if n < 0:
        raise ValueError("Number of points cannot be negative.")

    geod = Geodesic.WGS84

    line = geod.InverseLine(
        P1[0], P1[1],
        P2[0], P2[1]
    )

    total_distance = line.s13
    points = []

    for i in range(1, n + 1):
        fraction = i / (n + 1)
        distance = fraction * total_distance

        position = line.Position(
            distance,
            Geodesic.STANDARD
        )

        height = P1[2] + fraction * (P2[2] - P1[2])

        points.append({
            "latitude": position["lat2"],
            "longitude": position["lon2"],
            "height": height
        })

    return points