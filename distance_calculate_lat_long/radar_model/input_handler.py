def read_point(name):
    print(f"\nEnter {name} location:")

    latitude = float(input("Latitude (degrees): "))
    longitude = float(input("Longitude (degrees): "))
    height = float(input("Height (meters): "))

    if not -90 <= latitude <= 90:
        raise ValueError("Latitude must be between -90 and 90.")

    if not -180 <= longitude <= 180:
        raise ValueError("Longitude must be between -180 and 180.")

    return latitude, longitude, height