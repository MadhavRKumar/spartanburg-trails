import geopandas as gpd
import pandas as pd

import geoplot as gplt
import matplotlib.pyplot as plt

import folium as fl
from folium.plugins import Fullscreen, Geocoder

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

    road_length_minimum = 0.03  # minimum length of road in miles to consider for crash rate calculation
    collision_count_minimum = 1  # minimum number of collisions to consider for crash rate calculation

    if aadt > 0 and length_of_road_in_miles > road_length_minimum and collisions > collision_count_minimum:
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


def count_primary_factors(series):
    # sometimes we get one value, not a series, so we need to check for that
    if isinstance(series, str):
        series = pd.Series([series])

    if series.empty:
        series = pd.Series([])

    factors = series.value_counts().reset_index().rename(columns={0: 'count', 'index': 'PrimaryFactor'}).to_dict(orient='records')
    # get the factors with the highest count as formatted string
    # if there are multiple factors with the same count, return them all
    if not factors:
        return ''
    total_count = sum(factor['count'] for factor in factors)
    max_count = max(factor['count'] for factor in factors)

    # format the most_common_factors as a string with percentage
    formatted_factors = '<br/>'.join(f"{factor['PrimaryFactor']} ({(factor['count'] / total_count) * 100:.2f}%)" for factor in factors if factor['count'] == max_count)
    return formatted_factors if formatted_factors else ''


def calculate_city_crash_rates():
    collisions = read_file("spartanburg_collisions_2017-2021_clipped.geojson", "EPSG:32617")
    roads = read_file("complete_spartanburg_roads.geojson", "EPSG:32617")

    # separate collisions that have a JunctionType of 'Nonjunction' from all other collisions
    junction_collisions = collisions[~collisions['JunctionType'].isin(['Nonjunction', 'Unknown', 'Driveway'])]
    collisions = collisions[collisions['JunctionType'].isin(['Nonjunction', 'Unknown', 'Driveway'])]


    collisions['geometry'] = collisions['geometry'].buffer(10)
    

    # drop index_right 
    roads = roads.drop(columns=['index_right'], errors='ignore')

    collisions_per_road = roads.sjoin_nearest(collisions, distance_col="distance")

    def compare_street_names(row):
        if pd.isna(row['Street']) or pd.isna(row['STREET_NAM']):
            return 0.0
        return SequenceMatcher(None, row['Street'], row['STREET_NAM']).ratio()


    collisions_per_road['name_match'] = collisions_per_road.apply(compare_street_names, axis=1)
    # sort the collisions by name_match and drop duplicates
    collisions_per_road = collisions_per_road.sort_values(by='name_match', ascending=False).drop_duplicates(subset='AccidentNumber')



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
        collisions=('AccidentNumber', lambda x: ', '.join(x.unique())),
        primary_factors=('PrimaryFactor', count_primary_factors),
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
    # print(roads['length_of_road'].describe())
    ax = roads.plot(column='crash_rate_per_million_vmt', cmap='Wistia', legend=True, figsize=(20, 20)) 

    city_limits = read_file("spartanburg_city_limits.geojson")
    city_limits = city_limits.to_crs(4326)  # Set the CRS to WGS 84


    # print the top 10 roads with the highest crash rate
    # print(roads[['FactoredAA','BeginMileP', 'EndMilePoi', 'length_of_road', 
    #               'collision_count', 'RoadName', 'RouteLRS', 'crash_rate_per_million_vmt']].sort_values(by='crash_rate_per_million_vmt', ascending=False).head(10))


    # # save the roads with collision data to a GeoJSON file
    # # only save the columns we need
    roads_to_save = roads[['FactoredAA', 'BeginMileP', 'EndMilePoi', 'length_of_road', 'crash_rate_per_million_vmt', 
                           'collision_count', 'STREET_NAM', 'primary_factors', 'collisions',  'geometry']]

    roads_to_save.to_file("spartanburg_roads_with_collisions.geojson", driver="GeoJSON")


def calculate_crash_rates_for_junctions():
    collisions = read_file("spartanburg_collisions_2017-2021_clipped.geojson", "EPSG:32617")
    roads = read_file("complete_spartanburg_roads.geojson", "EPSG:32617")

    # remove collisions that have JunctionType in ['Nonjunction', 'Unknown', 'Driveway']
    junction_collisions = collisions[~collisions['JunctionType'].isin(['Nonjunction', 'Unknown', 'Driveway'])]

    # drop index_right 
    roads = roads.drop(columns=['index_right'], errors='ignore')


    # find all intersections of roads to find the junctions

    # buffer roads slightly to ensure that they intersect
    buffered_roads = roads.copy()
    buffered_roads['geometry'] = buffered_roads['geometry'].buffer(1)
    pre_dissolved_coords = buffered_roads.get_coordinates()
    dissolved_coords = buffered_roads.dissolve().get_coordinates()
    intersections = dissolved_coords.merge(pre_dissolved_coords,
                           on=["x","y"], how="left", indicator=True).query("_merge=='left_only'").drop_duplicates()

    intersections = gpd.GeoDataFrame(geometry=gpd.points_from_xy(x=intersections.x, y=intersections.y))

    # buffering causes the intersections to be 4 points, so we need to combine them into a single point
    intersections = intersections.dissolve().reset_index(drop=True)

    # convert geoms to list of point to use AgglomerativeClustering
    points = list(intersections['geometry'][0].geoms)

    clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=20, compute_full_tree=True)
    labels = clustering.fit_predict([[p.x, p.y] for p in points])

    # create a new GeoDataFrame with the clustered points
    clustered_points = []
    for label in set(labels):
        cluster_points = [points[i] for i in range(len(points)) if labels[i] == label]
        if cluster_points:
            clustered_points.append(MultiPoint(cluster_points).centroid)

    intersections = gpd.GeoDataFrame(geometry=clustered_points, crs=roads.crs)


    # join intersections with roads to calculate the AADT for each intersection
    buffered_roads['geometry'] = buffered_roads['geometry'].buffer(5)  # buffer roads slightly to ensure that they intersect
    intersections = intersections.sjoin_nearest(buffered_roads, distance_col="distance", how="left")
    # calculate the AADT for each intersection
    intersections['AADT'] = intersections['FactoredAA'].fillna(0)


    # combine the AADT for each intersection finding the sum of AADT and dividing by 2
    grouped_intersections = intersections.groupby(intersections.index).agg(
        intersection_AADT=('AADT', lambda x: x.sum() / 2),
        intersection_name=('RoadName', lambda x: ', '.join(x.unique())),
    ).reset_index(drop=True)

    intersections = intersections.merge(grouped_intersections, left_index=True, right_index=True, suffixes=('', '_mean'))

    intersections = intersections.drop(columns=['index_right'], errors='ignore')

    intersections = intersections.drop_duplicates(subset='geometry', keep='first').reset_index(drop=True)

    buffered_intersections = intersections.copy()
    buffered_intersections['geometry'] = buffered_intersections['geometry'].buffer(1)

    junctions_with_collisions = junction_collisions.sjoin_nearest(buffered_intersections , distance_col="distance", max_distance=25 )

    def get_unique_values_as_string(series):
        return ','.join(series.unique().astype(str)) if not series.empty else ''
    

    # aggregate the collisions by intersection and keep the AccidentNumber
    # only count unique AccidentNumbers
    junctions_with_collisions = junctions_with_collisions.groupby('index_right').agg(
        collision_count=('AccidentNumber', lambda x: x.nunique()),
        collisions=('AccidentNumber', lambda x: get_unique_values_as_string(x)),
        primary_factors=('PrimaryFactor', count_primary_factors),
    ).reset_index()

    # print(junctions_with_collisions.head())


    # merge the collision counts with the intersections GeoDataFrame
    intersections = intersections.merge(junctions_with_collisions, left_index=True, right_on='index_right', how='left')
    intersections['collision_count'] = intersections['collision_count'].fillna(0).astype(int)

    # drop any intersections with no collisions
    intersections = intersections[intersections['collision_count'] > 0]

    intersections['crash_rate_per_million_vmt'] = intersections['collision_count'] * 1e6 / (intersections['intersection_AADT'] * 365 * 5)


    # only keep the collision_count, AADT_mean, and geometry columns 
    intersections = intersections[['intersection_AADT', 'collision_count',  'crash_rate_per_million_vmt', 'collisions', 'primary_factors', 'intersection_name', 'geometry']].copy()

    # ax = roads.plot(figsize=(15,15), zorder=1, linewidth=1, color="gray")
    m = roads[['RoadName', 'FactoredAA', 'geometry']].to_crs(4326).explore(color='orange', name='Road')
    m = junction_collisions.to_crs(4326).explore(m=m, color='red', name='Collisions', marker_kwds={'radius': 2})
    m = intersections.to_crs(4326).explore(m=m, column='crash_rate_per_million_vmt', cmap='OrRd', name='Junctions with Collisions', marker_kwds={'radius': 5}, scheme='NaturalBreaks', legend=True)

    intersections.to_file("intersection_collisions.geojson", driver="GeoJSON")

def show_intersection_map():
    intersections = read_file("intersection_collisions.geojson", 4326)
    roads = read_file("complete_spartanburg_roads.geojson", 4326)

    m = fl.Map(max_bounds=True, zoom_start=12, min_zoom=12, max_zoom=30)
    m.fit_bounds([[34.89614214, -82.0035764], [34.98325887, -81.84675988]])
    roads.explore(m=m, color='orange', name='Roads', legend=True, figsize=(20, 20))
    intersections.explore(m=m, column='crash_rate_per_million_vmt', cmap='OrRd', name='Junctions with Collisions', marker_kwds={'radius': 5}, scheme='NaturalBreaks', legend=True)

    m.save("junctions_with_collisions.html")



def make_map():
    roads = read_file("spartanburg_roads_with_collisions.geojson", 4326)
    city_limits = read_file("spartanburg_city_limits.geojson", 4326)
    intersections = read_file("intersection_collisions.geojson", 4326)
    # collisions = read_file("spartanburg_collisions_2017-2021_clipped.geojson", 4326)

    # drop the collisions column from intersections
    intersections = intersections.drop(columns=['collisions'], errors='ignore')

    intersections = intersections.rename(columns={
        'intersection_AADT': 'Intersection AADT',
        'collision_count': 'Collisions',
        'crash_rate_per_million_vmt': 'Crash Rate',
        'intersection_name': 'Name',
        'primary_factors': 'Primary Factors',
    })

    # round the crash rate to 2 decimal places
    intersections['Crash Rate'] = intersections['Crash Rate'].round(2)
    intersections = intersections[['Name', 'Crash Rate','Intersection AADT', 'Collisions', 'Primary Factors', 'geometry']]


    # rename the columns to be more descriptive
    roads = roads.rename(columns={
        'FactoredAA': 'AADT',
        'BeginMileP': 'Begin Mile Point',
        'EndMilePoi': 'End Mile Point',
        'RoadName': 'Road Name',
        'length_of_road': 'Road Length (Miles)',
        'RouteLRS': 'Route LRS',
        'STREET_NAM': 'Name',
        'crash_rate_per_million_vmt': 'Crash Rate per Million VMT',
        'collision_count': 'Collisions',
        'primary_factors': 'Primary Factors',
    })
    roads['Crash Rate per Million VMT'] = roads['Crash Rate per Million VMT'].round(2)
    roads = roads[['Name', 'Crash Rate per Million VMT', 'Collisions', 'AADT', 'Road Length (Miles)', 'Primary Factors', 'geometry']]

    bounds = roads.total_bounds
    xmin, ymin, xmax, ymax = bounds

    m = fl.Map(min_lat=ymin, min_lon=xmin, max_lat=ymax, max_lon=xmax, zoom_start=12, min_zoom=12, max_zoom=30, zoom_control=False)

    fl.GeoJson(city_limits, style_function=lambda x: {'fillColor': 'black', 'color': 'black', 'weight': 2, 'fillOpacity': 0.4}, name='City Limits').add_to(m)

    m.fit_bounds([[ymin, xmin], [ymax, xmax]])

    roads.explore(m=m, column='Crash Rate per Million VMT', cmap='Reds', name='Crash Data', scheme='NaturalBreaks', legend=True, figsize=(20, 20),
                  style_kwds={'weight': 2})

    intersections.explore(m=m, column='Crash Rate', cmap='cool', name='Intersection Collisions', marker_kwds={'radius': 4}, scheme='NaturalBreaks', legend=True)

    # collisions.explore(m=m, style_kwds={'stroke': False, 'fillColor': 'black', 'fillOpacity': 0.5 }, name='Collisions', marker_kwds={'radius': 2})

    fl.LayerControl().add_to(m)

    m.add_child(Fullscreen(force_separate_button=True, position='topright'))
    m.add_child(Geocoder(position='topright', collapsed=True, zoom=12))

    m.save("spartanburg_roads_with_collisions.html")

    roads_to_print = roads[['Name', 'Crash Rate per Million VMT', 'Collisions', 'Primary Factors']]
    roads = roads[roads['Collisions'] > 5]

    intersections_to_print = intersections[['Name', 'Crash Rate', 'Collisions', 'Primary Factors']]
    intersections = intersections[intersections['Collisions'] > 5]

    # print(roads[['Name', 'Crash Rate per Million VMT', 'Collisions', 'Primary Factors']].sort_values(by='Crash Rate per Million VMT', ascending=False).head(10).to_markdown(index=False))
    #
    # print(intersections[['Name', 'Crash Rate', 'Collisions', 'Primary Factors']].sort_values(by='Crash Rate', ascending=False).head(10).to_markdown(index=False))



if __name__ == "__main__":
    # calculate_crash_rates_for_junctions()
    # calculate_city_crash_rates()
    make_map()

