import pymap3d as pm

# User se lat, lon, height lena (height meters mein)
def read_point(name):
    print(f"\n{name} ki location daalo:")
    lat = float(input("  Latitude  (deg): "))
    lon = float(input("  Longitude (deg): "))
    h = float(input("  Height    (m)  : "))
    return (lat, lon, h)


# Input
P1 = read_point("P1")
P2 = read_point("P2")
n = int(input("\nP1 aur P2 ke beech kitne points chahiye? "))

# 1 + 2. Slant Range, Azimuth, Elevation (P1 se naapa gaya)
az, el, R = pm.geodetic2aer(*P2, *P1)

print("\n--- Result ---")
print("Slant Range:", round(R, 2), "meters")
print("Azimuth:", round(az, 2), "degrees (north se)")
print("Elevation:", round(el, 2), "degrees (local horizon se)")

# 3. P1 aur P2 ke beech ke points (lat, lon, height)
print(f"\n{n} Intermediate Points:")

p1 = pm.geodetic2ecef(*P1)
p2 = pm.geodetic2ecef(*P2)

for i in range(1, n + 1):
    t = i / (n + 1)
    x, y, z = (a + t * (b - a) for a, b in zip(p1, p2))
    lat, lon, h = pm.ecef2geodetic(x, y, z)
    print(f"{i}. Lat: {lat:.6f}, Lon: {lon:.6f}, Height: {h:.2f} m")



# # P1: Radar
# lat1 = 28.6139
# lon1 = 77.2090
# h1 = 10.0

# # P2: Alag location par target
# lat2 = 28.6200
# lon2 = 77.2150
# h2 = 500.0