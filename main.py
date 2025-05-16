import os
from dotenv import load_dotenv
import geopandas as gpd

load_dotenv()


def read_file(path):
    try:
        gdf = gpd.read_file(path)
        gdf.geometry.set_crs("EPSG:4326", inplace=True)
        gdf["geometry"] = gdf["geometry"].to_crs("EPSG:32718")

        return gdf

    except Exception as e:
        print(f"Error reading file: {str(e)}")
        return None


def calculate_intersection_percent(gdf_polygon, gdf_polyline):
    # Get the total length of the polylines

    total_length = gdf_polyline.geometry.length.sum()

    # For each polygon, calculate what percentage of the total length is inside the polygon
    def calculate_percentage(polygon):
        # Create a mask for the lines that intersect with the polygon
        mask = gdf_polyline.intersects(polygon)
        # Calculate the total length of the lines that intersect with the polygon
        intersecting_length = gdf_polyline[mask].geometry.length.sum()
        # Calculate the percentage of the total length that is inside the polygon
        return (intersecting_length / total_length) * 100 if total_length > 0 else 0

    # Apply the function to each polygon and store the results in a new column
    gdf_polygon["percentage"] = gdf_polygon.geometry.apply(calculate_percentage)

    meters = gdf_polygon["percentage"] * total_length / 100
    miles = meters / 1609.344
    gdf_polygon["length"] = miles

    # Sort by District
    gdf_polygon = gdf_polygon.sort_values(by="District")

    # Return the updated GeoDataFrame with the new column
    return gdf_polygon[["District", "CouncilPer", "percentage", "length", "geometry"]]


if __name__ == "__main__":
    polygon_data = read_file(os.environ["POLYGON_DATA"])
    polyline_data = read_file(os.environ["POLYLINE_DATA"])

    output = calculate_intersection_percent(polygon_data, polyline_data)
    sum = output["percentage"].sum()

    output.geometry = output.geometry.to_crs("EPSG:4326")

    # output the results to a geojson file
    output.to_file("output.geojson", driver="GeoJSON")

    print(f"Total percentage of trails inside city districts: {sum:.2f}%")
