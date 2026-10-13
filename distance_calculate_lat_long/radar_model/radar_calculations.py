import pymap3d as pm


def calculate_radar_parameters(P1, P2):
    azimuth, elevation, slant_range = pm.geodetic2aer(
        P2[0], P2[1], P2[2],
        P1[0], P1[1], P1[2]
    )

    return {
        "slant_range": slant_range,
        "azimuth": azimuth,
        "elevation": elevation
    }