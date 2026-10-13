import numpy as np

# --- WGS84 Constants ---
A_SEMI_MAJOR = 6378137.0
FLATTENING = 1.0 / 298.257223563
B_SEMI_MINOR = A_SEMI_MAJOR * (1.0 - FLATTENING)
FIRST_ECCENTRICITY_SQ = FLATTENING * (2.0 - FLATTENING)


# --- Coordinate Conversions ---
def geodetic_to_ecef(lat, lon, h):
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)

    sin_lat = np.sin(lat_rad)
    cos_lat = np.cos(lat_rad)
    N = A_SEMI_MAJOR / np.sqrt(
        1.0 - FIRST_ECCENTRICITY_SQ * sin_lat**2
    )

    x = (N + h) * cos_lat * np.cos(lon_rad)
    y = (N + h) * cos_lat * np.sin(lon_rad)
    z = (N * (1.0 - FIRST_ECCENTRICITY_SQ) + h) * sin_lat

    return np.array([x, y, z], dtype=float)


def ecef_to_geodetic(x, y, z):
    p = np.sqrt(x**2 + y**2)

    # Handle the poles separately to avoid division by cos(latitude).
    if p < 1e-9:
        lat = 90.0 if z >= 0 else -90.0
        lon = 0.0
        h = abs(z) - B_SEMI_MINOR
        return lat, lon, h

    theta = np.arctan2(z * A_SEMI_MAJOR, p * B_SEMI_MINOR)

    lat_rad = np.arctan2(
        z + FIRST_ECCENTRICITY_SQ * B_SEMI_MINOR * np.sin(theta)**3,
        p - FIRST_ECCENTRICITY_SQ * A_SEMI_MAJOR * np.cos(theta)**3
    )
    lon_rad = np.arctan2(y, x)

    N = A_SEMI_MAJOR / np.sqrt(
        1.0 - FIRST_ECCENTRICITY_SQ * np.sin(lat_rad)**2
    )
    h = p / np.cos(lat_rad) - N

    return np.degrees(lat_rad), np.degrees(lon_rad), h


def local_enu_basis(lat, lon):
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)

    sin_lat, cos_lat = np.sin(lat_rad), np.cos(lat_rad)
    sin_lon, cos_lon = np.sin(lon_rad), np.cos(lon_rad)

    e_hat = np.array([-sin_lon, cos_lon, 0.0])
    n_hat = np.array([
        -sin_lat * cos_lon,
        -sin_lat * sin_lon,
        cos_lat
    ])
    u_hat = np.array([
        cos_lat * cos_lon,
        cos_lat * sin_lon,
        sin_lat
    ])

    return e_hat, n_hat, u_hat


# --- Case 1: Slant Range, Azimuth, Elevation ---
def compute_slant_range_az_el(p1_geo, p2_geo):
    if np.allclose(p1_geo, p2_geo):
        raise ValueError("ERR_COINCIDENT_POINTS")

    p1_ecef = geodetic_to_ecef(*p1_geo)
    p2_ecef = geodetic_to_ecef(*p2_geo)
    d_vector = p2_ecef - p1_ecef

    slant_range = np.linalg.norm(d_vector)
    e_hat, n_hat, u_hat = local_enu_basis(*p1_geo[:2])

    delta_e = np.dot(d_vector, e_hat)
    delta_n = np.dot(d_vector, n_hat)
    delta_u = np.dot(d_vector, u_hat)

    azimuth = np.degrees(np.arctan2(delta_e, delta_n)) % 360.0
    elevation = np.degrees(np.arctan2(
        delta_u, np.hypot(delta_e, delta_n)
    ))

    return slant_range, azimuth, elevation


def compute_intermediate_points(p1_geo, p2_geo, n_points=5):
    if n_points < 2:
        raise ValueError("n_points must be at least 2")

    p1_ecef = geodetic_to_ecef(*p1_geo)
    p2_ecef = geodetic_to_ecef(*p2_geo)

    points = []
    for i in range(n_points):
        s = i / (n_points - 1)
        p_s = p1_ecef + s * (p2_ecef - p1_ecef)
        points.append(ecef_to_geodetic(*p_s))

    return points


# --- Case 2: 3D Volume and Segment Intersection ---
def compute_segment_metrics(vertices_geo, r_geo, t_geo):
    if len(vertices_geo) != 8:
        raise ValueError("Provide 8 vertices for a 3D box.")

    # This case assumes the 8 vertices describe an axis-aligned
    # latitude/longitude/height box in geodetic coordinates.
    vertices = np.asarray(vertices_geo, dtype=float)
    lat_min, lat_max = vertices[:, 0].min(), vertices[:, 0].max()
    lon_min, lon_max = vertices[:, 1].min(), vertices[:, 1].max()
    h_min, h_max = vertices[:, 2].min(), vertices[:, 2].max()

    r_ecef = geodetic_to_ecef(*r_geo)
    t_ecef = geodetic_to_ecef(*t_geo)
    d_vec = t_ecef - r_ecef

    slant_range, azimuth, elevation = compute_slant_range_az_el(
        r_geo, t_geo
    )
    points = compute_intermediate_points(r_geo, t_geo, n_points=5)

    # Sample the straight ECEF segment to find portions inside the box.
    # This is a numerical approximation, not an exact curved-surface test.
    samples = 20001
    s_values = np.linspace(0.0, 1.0, samples)
    inside_s = []

    for s in s_values:
        pos = r_ecef + s * d_vec
        lat, lon, h = ecef_to_geodetic(*pos)

        if (lat_min <= lat <= lat_max and
                lon_min <= lon <= lon_max and
                h_min <= h <= h_max):
            inside_s.append(s)

    if inside_s:
        entry_s = min(inside_s)
        exit_s = max(inside_s)
        entry_geo = ecef_to_geodetic(*(r_ecef + entry_s * d_vec))
        exit_geo = ecef_to_geodetic(*(r_ecef + exit_s * d_vec))
        inside_length = (exit_s - entry_s) * slant_range
    else:
        entry_geo = None
        exit_geo = None
        inside_length = 0.0

    return {
        "A": r_geo,
        "B": t_geo,
        "R": slant_range,
        "Azimuth": azimuth,
        "Elevation": elevation,
        "Points": points,
        "Volume_Intersection": bool(inside_s),
        "Entry_Point": entry_geo,
        "Exit_Point": exit_geo,
        "Length_Inside_Volume_m": inside_length
    }


# --- Case 3: Angular Geometry and P3 Projection ---
def compute_case3_angular_geometry(p1_geo, p2_geo, p4_geo, d_vertical):
    _, az2, el2 = compute_slant_range_az_el(p1_geo, p2_geo)
    _, az4, el4 = compute_slant_range_az_el(p1_geo, p4_geo)

    p1_e = geodetic_to_ecef(*p1_geo)
    p2_e = geodetic_to_ecef(*p2_geo)
    p4_e = geodetic_to_ecef(*p4_geo)

    u2 = (p2_e - p1_e) / np.linalg.norm(p2_e - p1_e)
    u4 = (p4_e - p1_e) / np.linalg.norm(p4_e - p1_e)

    gamma = np.degrees(np.arccos(np.clip(np.dot(u2, u4), -1.0, 1.0)))
    delta_phi = (az4 - az2 + 180.0) % 360.0 - 180.0
    delta_theta = el4 - el2

    # Vertical offset is applied to ellipsoidal height, in metres.
    p3_geo = (p2_geo[0], p2_geo[1], p2_geo[2] - d_vertical)

    return {
        "Ray_P2": {"Azimuth": az2, "Elevation": el2},
        "Ray_P4": {"Azimuth": az4, "Elevation": el4},
        "Delta_Azimuth": delta_phi,
        "Delta_Elevation": delta_theta,
        "Included_Angle_Gamma": gamma,
        "P3": p3_geo
    }


# --- Main Execution ---
if __name__ == "__main__":
    print("=== GLV Module Execution Started ===")

    # Case 1
    p1 = (15.00, 73.00, 0.0)
    p2 = (15.10, 73.15, 5000.0)

    R, az, el = compute_slant_range_az_el(p1, p2)
    print(f"\n[Case 1] Slant Range: {R:.2f} m")
    print(f"Azimuth: {az:.4f} degrees")
    print(f"Elevation: {el:.4f} degrees")

    # Case 2: 3D box vertices
    vertices = [
        (15.02, 73.02, 0), (15.02, 73.10, 0),
        (15.08, 73.10, 0), (15.08, 73.02, 0),
        (15.02, 73.02, 3000), (15.02, 73.10, 3000),
        (15.08, 73.10, 3000), (15.08, 73.02, 3000)
    ]

    metrics = compute_segment_metrics(
        vertices,
        (14.98, 72.98, 0.0),
        (15.12, 73.14, 6000.0)
    )

    print(f"\n[Case 2] Segment Range: {metrics['R']:.2f} m")
    print(f"Azimuth: {metrics['Azimuth']:.4f} degrees")
    print(f"Elevation: {metrics['Elevation']:.4f} degrees")
    print(f"Intersects Volume: {metrics['Volume_Intersection']}")
    print(f"Entry Point: {metrics['Entry_Point']}")
    print(f"Exit Point: {metrics['Exit_Point']}")
    print(f"Length Inside Volume: {metrics['Length_Inside_Volume_m']:.2f} m")

    # Case 3
    case3 = compute_case3_angular_geometry(
        p1, p2, (15.15, 73.20, 6000.0), 1000.0
    )
    print(f"\n[Case 3] Angular Separation: "
          f"{case3['Included_Angle_Gamma']:.4f} degrees")
    print(f"P3: {case3['P3']}")

    print("\n=== Execution Completed Successfully ===")
