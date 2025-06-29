import os
from dotenv import load_dotenv
import geopandas as gpd

load_dotenv()


def read_file(path):
    try:
        gdf = gpd.read_file(path)
        gdf["geometry"] = gdf["geometry"].to_crs("EPSG:9822")

        return gdf

    except Exception as e:
        print(f"Error reading file: {str(e)}")
        return None





if __name__ == "__main__":
    cultural_district = read_file("https://umap.openstreetmap.fr/en/datalayer/1244274/9bee53a2-eb55-4a6f-8377-e7acfc1aeb32/?1751158845678")
    area_of_cultural_district = cultural_district.geometry.area.sum()

    non_places = read_file("https://umap.openstreetmap.fr/en/datalayer/1244274/15434ebd-91c7-4208-b2da-0eeffaa1e795/?1751114179678")
    places = read_file("https://umap.openstreetmap.fr/en/datalayer/1244274/4e21fc42-e669-4115-87a6-ede61c7a57bb/?1751158845677")

    non_places_in_cultural_district = non_places.clip(cultural_district)

    places_in_cultural_district = places.clip(cultural_district)
    places_in_cultural_district = places_in_cultural_district.overlay(non_places_in_cultural_district, how="difference")

    roads = cultural_district.overlay(places_in_cultural_district, how="difference")

    roads_and_non_places = roads.overlay(non_places_in_cultural_district, how="union")

    area_of_places_in_cultural_district = places_in_cultural_district.geometry.area.sum()

    print(f"Area of cultural district: {area_of_cultural_district:.2f} square meters")
    print(f"Area of places in cultural district: {area_of_places_in_cultural_district:.2f} square meters")

    area_of_non_places_in_cultural_district = roads_and_non_places.geometry.area.sum()
    print(f"Area of non-places in cultural district: {area_of_non_places_in_cultural_district:.2f} square meters")
    print()

    print(f"Ratio of places to non-places: {area_of_places_in_cultural_district / area_of_non_places_in_cultural_district:.2f}:1")

    print(f"Percentage of cultural district that is places: {area_of_places_in_cultural_district / area_of_cultural_district * 100:.2f}%")
    print(f"Percentage of cultural district that is non-places: {area_of_non_places_in_cultural_district / area_of_cultural_district * 100:.2f}%")

