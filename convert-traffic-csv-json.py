import csv

with open('traffic_count_data_2019.csv', 'r') as csv_file:
    csv_reader = csv.DictReader(csv_file)
    traffic_data = []
    
    for row in csv_reader:
        # create the LRS (Linear Referencing System) for each entry
        # Format: {county_id}{route_type}{route_number}{route_aux}
        # where county_id is 2 digits, route_type is 2 digits, route_number is 5 digits, route_aux is 2 digits
        
        # convert county, route designation, route number, and route auxiliary to integers
        row['County'] = int(row['County'])
        row['Route Designation'] = int(row['Route Designation'])
        row['Route Number'] = int(row['Route Number'])
        row['Route Auxiliary'] = int(row['Route Auxiliary'])

        lrs = f"{row['County']:02d}{row['Route Designation']:02d}{row['Route Number']:05d}{row['Route Auxiliary']:02d}"


        traffic_data.append({
            'lrs': lrs,
            'aadt': int(row['AADT']) if row['AADT'].isdigit() else 0,
            'location': row['Location'],
        })

    # Save the processed data to a JSON file
    import json
    with open('traffic_data_2019.json', 'w') as json_file:
        json.dump(traffic_data, json_file, indent=4)
