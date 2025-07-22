import geopandas as gpd
import pandas as pd

import geoplot as gplt
import matplotlib.pyplot as plt

import folium as fl
from folium.plugins import MeasureControl

from shapely import shortest_line, MultiPoint
from sklearn.cluster import AgglomerativeClustering

from difflib import SequenceMatcher

# non-stateroads - https://smpesri.scdot.org/arcgis/rest/services/EGIS_No_Imagery/MapServer/6/query
# statehighways - https://smpesi.scdot.org/arcgis/rest/services/EGIS_No_Imagery/MapServer/5/query
# state wide collision data (2017 - 2021) - https://services7.arcgis.com/TdsEnMqzMcTd7pnb/ArcGIS/rest/services/Aggregated_Collision_Points/FeatureServer/4/query


# case insensitive jq filter - https://siliconheaven.info/jq-filtering-by-properties/

# Steps:
# [x] Collect traffic count data for 2019 since it is the median year (as recommended by NCDOnT) 
# https://connect.ncdot.gov/resources/safety/Documents/TEAAS/Chapter%2008%20AADT.pdf
# [x] Extract aadt data for Spartanburg County 


# [x] build LRS for each entry in CSV and pair with AADT data
# 26 County ID, 07 Route Type, 00046 Route Number, 00 Route Aux, E Route Direction

# Pull roads (state and non-state) from SCDOT that are in Spartanburg 

# Import collision data for Spartanburg County from SCDOT (2017 - 2021)
# aggregate collisions by road segment

# For each road, calculate the crash rate per million vehicle miles traveled (VMT)
# using the AADT


def read_file(path, crs=3043):
    try:
        gdf = gpd.read_file(path)
        gdf["geometry"] = gdf["geometry"].to_crs(crs)  

        return gdf

    except Exception as e:
        print(f"Error reading file: {str(e)}")
        return None

def get_crash_data():
    collisions = None
    # query the server for the collision data
    # there are 20,000+ records in the dataset, so we will need to paginate through the results
    # by 1000 records at a time
    # colect them into one GeoDataFrame
    for i in range(0, 21):
        offset = i * 1000
        server_url = f"https://services7.arcgis.com/TdsEnMqzMcTd7pnb/ArcGIS/rest/services/Extrnl_Collisions/FeatureServer/16/query?where=County%3D%27Spartanburg%27+AND+RouteCategory%3C%3E%27Private+Property%27&objectIds=&geometry=-82.0035764%2C34.89614214%2C-81.84675988%2C34.98325887&geometryType=esriGeometryEnvelope&inSR=&spatialRel=esriSpatialRelIntersects&resultType=none&distance=0.0&units=esriSRUnit_Meter&relationParam=&returnGeodetic=false&outFields=*&returnGeometry=true&featureEncoding=esriDefault&multipatchOption=xyFootprint&maxAllowableOffset=&geometryPrecision=&outSR=&defaultSR=&datumTransformation=&applyVCSProjection=false&returnIdsOnly=false&returnUniqueIdsOnly=false&returnCountOnly=false&returnExtentOnly=false&returnQueryGeometry=false&returnDistinctValues=false&cacheHint=false&collation=&orderByFields=&groupByFieldsForStatistics=&outStatistics=&having=&resultOffset={offset}&resultRecordCount=&returnZ=false&returnM=false&returnTrueCurves=false&returnExceededLimitFeatures=true&quantizationParameters=&sqlFormat=none&f=pgeojson&token="

        response = gpd.read_file(server_url)
        if collisions is None:
            collisions = response
        else:
            collisions = pd.concat([collisions, response], ignore_index=True)

    # write the collisions to a GeoJSON file
    city_limits = read_file("spartanburg_city_limits.geojson")
    collisions = collisions.clip(city_limits)
    collisions.to_file("spartanburg_collisions_2017-2021_clipped.geojson", driver="GeoJSON")

def calculate_length_of_road(row):
    return float(row['EndMilePoi']) - float(row['BeginMileP'])


def calculate_crash_rate(row):
    aadt = row['FactoredAA']
    collisions = row['collision_count']
    length_of_road_in_miles = row['length_of_road']
    if aadt > 0 and length_of_road_in_miles > 0:
        return collisions * 1e6 / (aadt * 365 * 5 * length_of_road_in_miles)
    else:
        return 0

def get_city_and_county_roads():
    statewide_highways = read_file("https://services1.arcgis.com/VaY7cY9pvUYUP1Lf/ArcGIS/rest/services/Statewide_Highways/FeatureServer/0/query?where=1%3D1&objectIds=&geometry=-82.0035764%2C34.89614214%2C-81.84675988%2C34.98325887&geometryType=esriGeometryEnvelope&inSR=&spatialRel=esriSpatialRelIntersects&resultType=none&distance=0.0&units=esriSRUnit_Meter&relationParam=&returnGeodetic=false&outFields=*&returnGeometry=true&returnEnvelope=false&featureEncoding=esriDefault&multipatchOption=xyFootprint&maxAllowableOffset=&geometryPrecision=&outSR=&defaultSR=&datumTransformation=&applyVCSProjection=false&returnIdsOnly=false&returnUniqueIdsOnly=false&returnCountOnly=false&returnExtentOnly=false&returnQueryGeometry=false&returnDistinctValues=false&cacheHint=false&collation=&orderByFields=&groupByFieldsForStatistics=&outStatistics=&having=&resultOffset=&resultRecordCount=&returnZ=false&returnM=false&returnTrueCurves=false&returnExceededLimitFeatures=true&quantizationParameters=&sqlFormat=none&f=geojson&token=", "EPSG:32617")

    state_roads_with_aadt = read_file("https://smpesri.scdot.org/arcgis/rest/services/EGIS_No_Imagery/MapServer/6/query?where=1%3D1&text=&objectIds=&time=&geometry=-9128596.368%2C4149776.218%2C-9111139.6325%2C4161606.3241&geometryType=esriGeometryEnvelope&inSR=&spatialRel=esriSpatialRelIntersects&relationParam=&outFields=*&returnGeometry=true&returnTrueCurves=false&maxAllowableOffset=&geometryPrecision=&outSR=&having=&returnIdsOnly=false&returnCountOnly=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&returnZ=false&returnM=false&gdbVersion=&historicMoment=&returnDistinctValues=false&resultOffset=&resultRecordCount=&queryByDistance=&returnExtentOnly=false&datumTransformation=&parameterValues=&rangeValues=&quantizationParameters=&featureEncoding=esriDefault&f=geojson", "EPSG:32617")

    spartanburg_city_roads = read_file("https://smpesri.scdot.org/arcgis/rest/services/EGIS_No_Imagery/MapServer/5/query?where=1%3D1&text=&objectIds=&time=&geometry=-9128596.368%2C4149776.218%2C-9111139.6325%2C4161606.3241&geometryType=esriGeometryEnvelope&inSR=&spatialRel=esriSpatialRelIntersects&relationParam=&outFields=*&returnGeometry=true&returnTrueCurves=false&maxAllowableOffset=&geometryPrecision=&outSR=&having=&returnIdsOnly=false&returnCountOnly=false&orderByFields=&groupByFieldsForStatistics=&outStatistics=&returnZ=false&returnM=false&gdbVersion=&historicMoment=&returnDistinctValues=false&resultOffset=&resultRecordCount=&queryByDistance=&returnExtentOnly=false&datumTransformation=&parameterValues=&rangeValues=&quantizationParameters=&featureEncoding=esriDefault&f=geojson","EPSG:32617")


    # buffer the geometry of the roads to ensure that they are not too small
    state_roads_with_aadt['geometry'] = state_roads_with_aadt['geometry'].buffer(15)

    # some roads don't have AADT data so we spataially join the roads with AADT data
    roads = statewide_highways.sjoin_nearest(state_roads_with_aadt, distance_col="distance", how="left", max_distance=50)

    # for any duplicated indices keep the ones where RoadName == STREET_NAM
    roads['name_match'] = roads.apply(lambda row: SequenceMatcher(None, row['RoadName'], row['STREET_NAM']).ratio(), axis=1)
    roads = roads.sort_values(by='name_match', ascending=False).drop_duplicates(subset='FID')

    # roads = roads.drop(columns=[col for col in roads.columns if col.endswith('_right')], errors='ignore')

    # rename the columns to remove the _left suffix
    roads = roads.rename(columns=lambda x: x.replace('_left', ''))

    #drop index_right

    # merge the BEG_MILEPO and BeginMimerge columns after joinleP columns into BeginMileP
    roads['BeginMileP'] = roads['BEG_MILEPO'].fillna(roads['BeginMileP'])
    roads = roads.drop(columns=['BEG_MILEPO'], errors='ignore')

    # merge the END_MILEPO and EndMilePoi columns into EndMilePoi
    roads['EndMilePoi'] = roads['END_MILEPO'].fillna(roads['EndMilePoi'])
    roads = roads.drop(columns=['END_MILEPO'], errors='ignore')


    roads['length_of_road'] = roads.apply(calculate_length_of_road, axis=1)

    roads = roads.merge(state_roads_with_aadt[['geometry']], left_on='index_right', right_index=True, suffixes=('', '_right'), how='left')

    connectors = roads.apply(lambda row: shortest_line(row['geometry'], row['geometry_right']), axis=1)
    connectors = gpd.GeoDataFrame(geometry=connectors, crs=roads.crs)
    connectors['distance'] = connectors.geometry.length

    roads = roads.drop(columns=['geometry_right', 'index_right', 'name_match'], errors='ignore')

    # combine the state highways and the city roads
    roads = pd.concat([roads, spartanburg_city_roads], ignore_index=True)

    # fill STREET_NAM with RoadName if it is NaN
    roads['STREET_NAM'] = roads['STREET_NAM'].fillna(roads['RoadName'])

    city_limits = read_file("spartanburg_city_limits.geojson", "EPSG:32617")

    roads = roads.clip(city_limits)

    # clip the roads to the city limits
    # roads = roads.clip(city_limits)
    roads = roads.to_crs(4326)  # Set the CRS to WGS 84Roads


    m = statewide_highways.to_crs(4326).explore(color='blue', name='Statewide Highways')
    m = state_roads_with_aadt.to_crs(4326).explore(m=m, color='red', name='State Roads with AADT') 
    # m = spartanburg_city_roads.to_crs(4326).explore(m=m, color='green', name='Spartanburg City Roads')
    m = connectors.to_crs(4326).explore(m=m, color='purpfrom difflib import SequenceMatcherle', name='Connectors') 
    m = roads.explore(m=m, color='orange', name='All Spartanburg Roads with AADT')

    fl.LayerControl().add_to(m)
    m.save("spartanburg_roads.html")

    # save the roads to a GeoJSON file
    roads.to_file("complete_spartanburg_roads.geojson", driver="GeoJSON")

def calculate_city_crash_rates():
    collisions = read_file("spartanburg_collisions_2017-2021_clipped.geojson", "EPSG:32617")
    roads = read_file("complete_spartanburg_roads.geojson", "EPSG:32617")

    # separate collisions that have a JunctionType of 'Nonjunction' from all other collisions
    junction_collisions = collisions[collisions['JunctionType'] == 'Nonjunction']
    collisions = collisions[collisions['JunctionType'] != 'Nonjunction']


    collisions['geometry'] = collisions['geometry'].buffer(1)
    

    # drop index_right 
    roads = roads.drop(columns=['index_right'], errors='ignore')

    collisions_per_road = roads.sjoin_nearest(collisions, distance_col="distance")

    # function to find most common street name from agg, if there are multiple names
    # with the same count, it will return the first one
    def most_common_street_name(street_names):
        if street_names.empty:
            return None
        # count the occurrences of each street name
        counts = street_names.value_counts()
        # return the most common street name
        return counts.idxmax() if not counts.empty else None

    def get_unique_values(series):
        return series.unique().tolist() if not series.empty else []

    # aggregate the collisions by road segment and keep the Street Name and OBJECTID
    # store all of the collision data for each road segment
    collisions_per_road = collisions_per_road.groupby(collisions_per_road.index).agg(
        collision_count=('AccidentNumber', 'size'),
        collisions=('AccidentNumber', lambda x: list(x)),
    ).reset_index()


    # merge the collision counts with the roads GeoDataFrame
    roads = roads.merge(collisions_per_road, left_index=True, right_on='index', how='left')


    roads['collision_count'] = roads['collision_count'].fillna(0).astype(int)
    roads['length_of_road'] = roads.apply(calculate_length_of_road, axis=1)
    roads['crash_rate_per_million_vmt'] = roads.apply(calculate_crash_rate, axis=1)

    # only keep roads with a crash rate greater than 1
    # roads = roads[roads['collision_count'] > 1]

    # drop the first quarter of the roads with the lowest crash rate


    roads['geometry'] = roads['geometry'].to_crs(4326)  # Set the CRS to WGS 84
    print(roads['length_of_road'].describe())
    roads = roads[roads['length_of_road'] > 0.03]

    ax = roads.plot(column='crash_rate_per_million_vmt', cmap='Wistia', legend=True, figsize=(20, 20)) 

    city_limits = read_file("spartanburg_city_limits.geojson")
    city_limits = city_limits.to_crs(4326)  # Set the CRS to WGS 84
    gplt.polyplot(city_limits.explode(), edgecolor='black', facecolor='black', linewidth=1, ax=ax)

    plt.show()


    # print the top 10 roads with the highest crash rate
    print(roads[['FactoredAA','BeginMileP', 'EndMilePoi', 'length_of_road', 
                  'collision_count', 'RoadName', 'RouteLRS', 'crash_rate_per_million_vmt']].sort_values(by='crash_rate_per_million_vmt', ascending=False).head(10))


    # # save the roads with collision data to a GeoJSON file
    # # only save the columns we need
    roads_to_save = roads[['FactoredAA', 'BeginMileP', 'EndMilePoi', 'length_of_road', 'crash_rate_per_million_vmt', 
                           'collision_count', 'STREET_NAM',  'geometry']]
    roads_to_save.to_file("spartanburg_roads_with_collisions.geojson", driver="GeoJSON")


def calculate_crash_rates_for_junctions():
    collisions = read_file("spartanburg_collisions_2017-2021_clipped.geojson", "EPSG:32617")
    roads = read_file("complete_spartanburg_roads.geojson", "EPSG:32617")

    # remove collisions that have JunctioType in ['Nonjunction', 'Unknown', 'Driveway']
    junction_collisions = collisions[~collisions['JunctionType'].isin(['Nonjunction', 'Unknown', 'Driveway'])]

    # drop index_right 
    roads = roads.drop(columns=['index_right'], errors='ignore')


    # find all intersections of roads to find the junctions

    # buffer roads slightly to ensure that they intersect
    buffered_roads = roads.copy()
    buffered_roads['geometry'] = buffered_roads['geometry'].buffer(0.5)
    pre_dissolved_coords = buffered_roads.get_coordinates()
    dissolved_coords = buffered_roads.dissolve().get_coordinates()
    intersections = dissolved_coords.merge(pre_dissolved_coords,
                           on=["x","y"], how="left", indicator=True).query("_merge=='left_only'").drop_duplicates()

    intersections = gpd.GeoDataFrame(geometry=gpd.points_from_xy(x=intersections.x, y=intersections.y))

    # buffering causes the intersections to be 4 points, so we need to combine them into a single point
    intersections = intersections.dissolve().reset_index(drop=True)

    # convert geoms to list of point to use AgglomerativeClustering
    points = list(intersections['geometry'][0].geoms)

    clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=30, compute_full_tree=True)
    labels = clustering.fit_predict([[p.x, p.y] for p in points])

    # create a new GeoDataFrame with the clustered points
    clustered_points = []
    for label in set(labels):
        cluster_points = [points[i] for i in range(len(points)) if labels[i] == label]
        if cluster_points:
            clustered_points.append(MultiPoint(cluster_points).centroid)

    intersections = gpd.GeoDataFrame(geometry=clustered_points, crs=roads.crs)


    # join intersections with roads to calculate the AADT for each intersection
    intersections = intersections.sjoin_nearest(buffered_roads, distance_col="distance", how="left")
    # calculate the AADT for each intersection
    intersections['AADT'] = intersections['FactoredAA'].fillna(0)


    # combine the AADT for each intersection finding the sum of AADT and dividing by 2
    grouped_intersections = intersections.groupby(intersections.index).agg(
        intersection_AADT=('AADT', lambda x: x.sum() / 2),
    ).reset_index(drop=True)

    intersections = intersections.merge(grouped_intersections, left_index=True, right_index=True, suffixes=('', '_mean'))
    intersections = intersections.drop(columns=['index_right'], errors='ignore')


    buffered_intersections = intersections.copy()
    buffered_intersections['geometry'] = buffered_intersections['geometry'].buffer(1)

    junctions_with_collisions = junction_collisions.sjoin_nearest(buffered_intersections , distance_col="distance")

    # aggregate the collisions by intersection and keep the AccidentNumber
    junctions_with_collisions = junctions_with_collisions.groupby('index_right').agg(
        collision_count=('AccidentNumber', 'size'),
        # collisions=('AccidentNumber', lambda x: list(x)),
    ).reset_index()

    # merge the collision counts with the intersections GeoDataFrame
    intersections = intersections.merge(junctions_with_collisions, left_index=True, right_on='index_right', how='left')

    # only keep the collision_count, AADT_mean, and geometry columns 
    intersections = intersections[['intersection_AADT', 'collision_count', 'geometry']].copy()


    # ax = roads.plot(figsize=(15,15), zorder=1, linewidth=1, color="gray")
    m = roads.to_crs(4326).explore(color='orange', name='Road')
    m = junction_collisions.to_crs(4326).explore(m=m, color='red', name='Collisions', marker_kwds={'radius': 2})
    m = intersections.to_crs(4326).explore(m=m, column='collision_count', name='Junctions with Collisions', marker_kwds={'radius': 5})
    m.add_child(MeasureControl())
    m.save("junctions_with_collisions.html")

    # junctions_with_collisions.plot(ax=ax, zorder=3, color="blue", markersize=5)


def make_map():
    roads = read_file("spartanburg_roads_with_collisions.geojson", 4326)

    # rename the columns to be more descriptive
    roads = roads.rename(columns={
        'FactoredAA': 'AADT',
        'BeginMileP': 'Begin Mile Point',
        'EndMilePoi': 'End Mile Point',
        'RoadName': 'Road Name',
        'RouteLRS': 'Route LRS',
        'crash_rate_per_million_vmt': 'Crash Rate per Million VMT',
        'collision_count': 'Collision Count',
    })

    bounds = roads.total_bounds
    xmin, ymin, xmax, ymax = bounds
    m = fl.Map(max_bounds=True, min_lat=ymin, min_lon=xmin, max_lat=ymax, max_lon=xmax, zoom_start=12, min_zoom=12, max_zoom=18)
    m.fit_bounds([[ymin, xmin], [ymax, xmax]])
    roads.explore(m=m,column='Crash Rate per Million VMT', cmap='YlOrRd', name='Crash Data', legend=True, figsize=(20, 20))
    fl.TileLayer('CartoDB.DarkMatter').add_to(m)
    fl.LayerControl().add_to(m)
    m.save("spartanburg_roads_with_collisions.html")



if __name__ == "__main__":
    calculate_crash_rates_for_junctions()
