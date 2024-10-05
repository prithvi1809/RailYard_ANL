import os
import sys
import copy
import yaml
import time
import random
import requests
import osm2gmns as og
import osm2rail as orl
import subprocess, re
from DLSim import DLSim
from .painter import painter
from grid2demand import GRID2DEMAND
import datetime
import matplotlib
import folium
from folium.plugins import FloatImage
import urllib.request
from PIL import Image
import math
import matplotlib.pyplot as plt
# map html-style visualization
import urllib.request
from PIL import Image

from shapely.wkt import loads as load_wkt
import uuid

import json
import pandas as pd
import numpy as np
import shapely.geometry as geometry
from shapely import wkt, Point, LineString
import ast
from shapely.geometry import MultiLineString, LineString, Polygon
from shapely.ops import unary_union

matplotlib.use('TkAgg')

config_path = "/home/local/ASURITE/longchao/Desktop/project/LLM4Traffic/OpenTI/pivotalAgent/Configs/path.yaml"
api_config = "/home/local/ASURITE/longchao/Desktop/project/LLM4Traffic/OpenTI/pivotalAgent/Configs/config.yaml"
with open(config_path, 'r') as file:
    tsc_root = yaml.safe_load(file)["LibSignal"]["tsc_root"]

with open(config_path, 'r') as file:
    filter_store_base = yaml.safe_load(file)["FIlter_Store_Base"]

def func_prompt(name, description):
    def decorator(func):
        func.name = name
        func.description = description
        return func

    return decorator

class ask4Area:
    def __init__(self) -> None:
        pass

    @func_prompt(name="queryAreaRange",
             description="""
             This tool is used to obtain the area information of a interest point on map.
             Consider using this tool when asked "Where's the interest point?"
             The output will tell you whether you have finished this command successfully. 
             """)

    def embody(self, target: str) -> str:

        def get_bounding_box(query):
            # URL for the API request with the provided query
            url = f"https://nominatim.openstreetmap.org/search.php?q={query}&polygon_geojson=1&format=jsonv2"

            # Headers to simulate the request
            headers = {
                'Accept': '*/*',
                'Accept-Encoding': 'gzip, deflate, br, zstd',
                'Accept-Language': 'en-US,en;q=0.9',
                'Referer': 'https://nominatim.openstreetmap.org/ui/search.html?q=arizona+state+university+',
                'Sec-Ch-Ua': '"Google Chrome";v="119", "Chromium";v="119", "Not?A_Brand";v="24"',
                'Sec-Ch-Ua-Mobile': '?0',
                'Sec-Ch-Ua-Platform': '"Linux"',
                'Sec-Fetch-Dest': 'empty',
                'Sec-Fetch-Mode': 'cors',
                'Sec-Fetch-Site': 'same-origin',
                'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/119.0.0.0 Safari/537.36'
            }

            try:
                # Make the GET request with headers
                response = requests.get(url, headers=headers)

                if response.status_code == 200:
                    data = response.json()

                    # Assuming we need the bounding box information from the first result
                    if len(data) > 0:
                        bounding_box = data[0].get('boundingbox', [])
                        if bounding_box:
                            return bounding_box
                        else:
                            print("Bounding box not found in the response.")
                    else:
                        print("No results found in the response.")
                else:
                    print(f"Failed to get a response. Status code: {response.status_code}")
            except Exception as e:
                print(f"An error occurred: {e}")

            return None
        
        location = get_bounding_box(target)

        try:
            return f"If you successuflly loacated the map. just return You have successfully located the map of: {target}. And your final answer should include this sentence without changing anything: The longitude and latitude area of interested {target} is:{location}."
        except json.JSONDecodeError:
            return "Invalid location format. Please provide a valid list of longitude and latitude values."

class GoogleMapsLayers:
        ROADMAP = "v"
        TERRAIN = "p"
        ALTERED_ROADMAP = "r"
        SATELLITE = "s"
        TERRAIN_ONLY = "t"
        HYBRID = "y"

class GoogleMapDownloader:
        def __init__(self, min_long, min_lat, max_long, max_lat, zoom=12, layer=GoogleMapsLayers.ALTERED_ROADMAP):
            self._west = min_long
            self._east = max_long
            self._south = min_lat
            self._north = max_lat
            self._zoom = zoom
            self._layer = layer

        def getXY(self):
            tile_size = 256
            num_tiles = 1 << self._zoom

            top_left_x = (tile_size / 2 + self._west * tile_size / 360.0) * num_tiles // tile_size
            sin_top_left_y = math.sin(self._north * (math.pi / 180.0))
            top_left_y = ((tile_size / 2) + 0.5 * math.log((1 + sin_top_left_y) / (1 - sin_top_left_y)) * -(
                tile_size / (2 * math.pi))) * num_tiles // tile_size

            bottom_right_x = (tile_size / 2 + self._east * tile_size / 360.0) * num_tiles // tile_size
            sin_bottom_right_y = math.sin(self._south * (math.pi / 180.0))
            bottom_right_y = ((tile_size / 2) + 0.5 * math.log((1 + sin_bottom_right_y) / (1 - sin_bottom_right_y)) * -(
                tile_size / (2 * math.pi))) * num_tiles // tile_size

            return int(top_left_x), int(top_left_y), int(bottom_right_x), int(bottom_right_y)

        def generateImage(self, **kwargs):
            start_x = kwargs.get('start_x', None)
            start_y = kwargs.get('start_y', None)
            tile_width = kwargs.get('tile_width', 5)
            tile_height = kwargs.get('tile_height', 5)

            if start_x is None or start_y is None:
                start_x, start_y, _, _ = self.getXY()

            width, height = 256 * tile_width, 256 * tile_height
            map_img = Image.new('RGB', (width, height))

            for x in range(0, tile_width):
                for y in range(0, tile_height):
                    url = f'https://mt0.google.com/vt?lyrs={self._layer}&x=' + str(start_x + x) + '&y=' + str(start_y + y) + '&z=' + str(
                        self._zoom)

                    current_tile = str(x) + '-' + str(y)
                    urllib.request.urlretrieve(url, current_tile)

                    im = Image.open(current_tile)
                    map_img.paste(im, (x * 256, y * 256))
                    os.remove(current_tile)

            return map_img 

class GoogleMapsLayers:
        ROADMAP = "v"
        TERRAIN = "p"
        ALTERED_ROADMAP = "r"
        SATELLITE = "s"
        TERRAIN_ONLY = "t"
        HYBRID = "y"

class GoogleMapDownloader:
        def __init__(self, min_longtitude, min_latitude, max_longtitude, max_latitude, zoom=12, layer=GoogleMapsLayers.ALTERED_ROADMAP):
            self._west = min_longtitude
            self._east = max_longtitude
            self._south = min_latitude
            self._north = max_latitude
            self._zoom = zoom
            self._layer = layer

        def getXY(self):
            tile_size = 256
            num_tiles = 1 << self._zoom

            top_left_x = (tile_size / 2 + self._west * tile_size / 360.0) * num_tiles // tile_size
            sin_top_left_y = math.sin(self._north * (math.pi / 180.0))
            top_left_y = ((tile_size / 2) + 0.5 * math.log((1 + sin_top_left_y) / (1 - sin_top_left_y)) * -(
                tile_size / (2 * math.pi))) * num_tiles // tile_size

            bottom_right_x = (tile_size / 2 + self._east * tile_size / 360.0) * num_tiles // tile_size
            sin_bottom_right_y = math.sin(self._south * (math.pi / 180.0))
            bottom_right_y = ((tile_size / 2) + 0.5 * math.log((1 + sin_bottom_right_y) / (1 - sin_bottom_right_y)) * -(
                tile_size / (2 * math.pi))) * num_tiles // tile_size

            return int(top_left_x), int(top_left_y), int(bottom_right_x), int(bottom_right_y)

        def generateImage(self, **kwargs):
            start_x = kwargs.get('start_x', None)
            start_y = kwargs.get('start_y', None)
            tile_width = kwargs.get('tile_width', 5)
            tile_height = kwargs.get('tile_height', 5)

            if start_x is None or start_y is None:
                start_x, start_y, _, _ = self.getXY()

            width, height = 256 * tile_width, 256 * tile_height
            map_img = Image.new('RGB', (width, height))

            for x in range(0, tile_width):
                for y in range(0, tile_height):
                    url = f'https://mt0.google.com/vt?lyrs={self._layer}&x=' + str(start_x + x) + '&y=' + str(start_y + y) + '&z=' + str(
                        self._zoom)

                    current_tile = str(x) + '-' + str(y)
                    urllib.request.urlretrieve(url, current_tile)

                    im = Image.open(current_tile)
                    map_img.paste(im, (x * 256, y * 256))

                    os.remove(current_tile)

            return map_img 

class showMap:
    def __init__(self) -> None:
        pass 
        # Specify the bounding box coordinates

    @func_prompt(name="showOnMap",
            description="""
             This tool is used to show the interested loaction on Map, like the ASU campus area on map.
             Consider using this tool when asked "Can you show me this area on the city map ?"
             The output will tell you whether you have finished this command successfully.
             """)
    def embody(self, target: str) -> str:
        
        with open(api_config, 'r') as file:
            google_api = yaml.safe_load(file)["Google_Map"]
        api_key = google_api
        target = eval(target)
        for i in range(len(target)):
            target[i] = float(target[i])

        min_lat, max_lat, min_long, max_long = target[0],target[1],target[2],target[3]
        center_lat = (min_lat + max_lat) / 2
        center_long = (min_long + max_long) / 2
        center = f'{center_lat},{center_long}'

        zoom = 13  # You might need to adjust this value based on the area covered by your bounding box

        # Define the URL for the Static Map API
        url = f'https://maps.googleapis.com/maps/api/staticmap?center={center}&zoom={zoom}&size=600x400&key={api_key}'

        # Make the request and get the response
        response = requests.get(url)

        def create_and_show_map(bounding_box):
            min_lat, max_lat, min_lon, max_lon = bounding_box[0], bounding_box[1],bounding_box[2],bounding_box[3]
            # Center of the map
            center_lat = (min_lat + max_lat) / 2
            center_lon = (min_lon + max_lon) / 2

            # Create a map object
            map_obj = folium.Map(location=[center_lat, center_lon], zoom_start=13)

            # Add bounding box as a rectangle
            folium.Rectangle(
                bounds=[[min_lat, min_lon], [max_lat, max_lon]],
                color='#ff7800',
                fill=True,
                fill_color='#ffff00',
                fill_opacity=0.2
            ).add_to(map_obj)

            return map_obj
        with open(config_path, 'r') as file:
            data = yaml.safe_load(file)
        target_map_png = data['showMap']['targetMapPng']
        html_path = data['showMap']['htmlPath']
        try:
            response = requests.get(url)
            if response.status_code == 200:
                with open(target_map_png, 'wb') as file:
                    file.write(response.content)
            else:
                print(f"Failed to get the map. Status code: {response.status_code}")
        except:
            return "Error for finding the place"

        
        map = create_and_show_map(target)
        map.save(html_path)

        return f"You have successfully find the map of: {target}. And your final answer should include this sentence without changing anything: The map area of interested {target} is: `{target_map_png}`, and dynamic version please check the link above:`{html_path}`."

class autoDownloadNetwork:
    def __init__(self, base_loc: str) -> None:
        self.base_loc = base_loc
    
    @func_prompt(name="autoDownloadOpenStreetMapFile",
             description="""
             This tool is used to automatically download the certain area's map data from OpenStreetMap. The downloaded file is data ends with '.osm'.
             This tool will return the file path that the downloaded file has been stored. 
             Consider using this tool if the question is anything about get me the map of area, or download the map from openstreetmap, or I want to download the map data from openstreetmap(OpenStreetMap).
             The output will tell you whether you have finished this command successfully.
             """
             )

    def embody(self, target: str) -> str:

        # try: 
            print("get the data in download: {}".format(target))
            desired = target.replace(" ", "").split(",")[-1].replace('"', '').replace("[", "").replace("]", "")
            print("desired:"+desired)
            # min_long, min_lat, max_long, max_lat = target.replace("[", "").replace("]", "").strip().replace(" ", "").split(",")[:4]
            coordinates, _ = target.split(',')

            # Then, split the coordinates by '*' to get the individual values
            min_lat, max_lat, min_long, max_long = coordinates.split('*')

            print("long_min, lat_min, long_max, lat_max")
            print(min_lat, max_lat, min_long, max_long)
            url = "https://www.openstreetmap.org/api/0.6/map?bbox={}%2C{}%2C{}%2C{}".format(min_long, min_lat, max_long, max_lat)
            print("url:")
            # min_lon,min_lat,max_lon,max_lat.
            print(url)
            response = requests.get(url)

            with open(config_path, 'r') as file:
                osm_base = yaml.safe_load(file)["OSM_BASE"]
            if response.status_code == 200:
                file_path = osm_base + desired
                # file_path = self.base_loc+ desired
                with open(file_path, 'wb') as file:
                    file.write(response.content)
                if not os.path.exists(file_path):
                    return f"The requested cannot be successfully downloaded."
            else:
                file_path = 'None'
                print("Failed to retrieve the data.")
                return f"The requested cannot be successfully downloaded."
           
            time.sleep(2)
            return f"The requested have successfully downloaded and saved at: {file_path}. And your final answer should include this sentence without changing anything: The file saved location is at: `{file_path}`."
        # except FileNotFoundError as e:
        #     return f"The requested cannot be successfully downloaded because your request was too large. Either request a smaller area, or use planet.osm."

def extract_filepath(s):
    match = re.search(r'data/.*\.log', s)
    return match.group(0) if match else None

def update_episode_in_base_yaml(episode_value: str) -> None:
    with open(config_path, 'r') as file:
        data = yaml.safe_load(file)
    base_yaml_path = data['LibSignal']['tsc_base']

    with open(base_yaml_path, 'r') as file:
        data = yaml.safe_load(file)

    # Update the episodes value
    data["trainer"]["episodes"] = int(episode_value)

    with open(base_yaml_path, 'w') as file:
        yaml.dump(data, file, default_flow_style=False)

class simulateOnLibSignal:
    def __init__(self, base_dir: str, save_dir: str) -> None:
        self.base_dir = base_dir
        self.save_dir = save_dir
        self.temp_directfeedback = tsc_root
    
    @func_prompt(name="simulateOnLibsignal",
             description="""
             This tool is used to execute the simulation on open source library called LibSignal. 
             Consider using this tool when it is mentioned of LibSignal and run/execution/simulate in/on LibSignal.
             The output will tell you whether you have finished this command successfully.
             """)
    
    
    def embody(self, target: str) -> str:
        print("target:")
        print(target)

        support_simulators = ["cityflow", "sumo"]
        support_algorithms_1 = ["dqn", "frap", "presslight", "mplight"]
        support_algorithms_2 = ["fixedtime", "sotl"]
        # try:
        print("get the command:{}".format(target))


        # Using regular expressions to extract information
        matches = re.match(r'([^,]+),\s*([^,]+),\s*(\d+)', target)

        if matches:
            simulator = matches.group(1).strip()
            algorithm = matches.group(2).strip()
            episode = int(matches.group(3).strip())
            print(simulator, algorithm, episode)
        else:
            # Using string manipulation to extract information
            simulator = None
            algorithm = None
            episode = None

            for param in target.replace(" ", "").split(","):
                if param.startswith("simulator="):
                    simulator = param.split("=")[1]
                elif param.startswith("algorithm="):
                    algorithm = param.split("=")[1]
                elif param.startswith("episode="):
                    episode = int(param.split("=")[1])

        print(simulator, algorithm, episode)
        
       
            
        # Set default value if episode is not provided
        episode = int(episode) if episode is not None else 5
        
        # Call the function to update the episode in base.yml
        update_episode_in_base_yaml(episode)

        
        with open(config_path, 'r') as file:
            log_example = yaml.safe_load(file)["LibSignal"]["log_example"]
            
        with open(config_path, 'r') as file2:
            libsignal_root = yaml.safe_load(file2)["LibSignal"]["root_path"]

        #successful run to here
        print(simulator, algorithm, episode)
        if simulator.lower() == "cityflow":
            if algorithm.lower() in support_algorithms_1:

           
                base_path = libsignal_root
                print(base_path)
                
                # Generate a timestamp for the log file
                timestamp = datetime.datetime.now().strftime("%Y_%m_%d-%H_%M_%S")
                print(timestamp)
                # Define the relative path for the log file
                log_file_relative_path = "data/output_data/tsc/cityflow_{}/cityflow1x1/0/logger/{}_DTL.log".format(algorithm,timestamp)
                print(log_file_relative_path)
                # # Combine the base path and relative path to get the full log file path
                log_file_path = os.path.join(base_path, log_file_relative_path)
                print(log_file_path)
                completed_process = subprocess.run(["python3", base_path+"run.py", "-a", algorithm,"-tt", timestamp], capture_output=True, text=True, cwd=base_path)
                # Check if the subprocess completed successfully
                output = completed_process.stdout
                print("output:")
                print(output)
                
                read_path = log_file_path
                print(read_path)

                saved_image = painter({'hz1x1': read_path}, ['epoch', 'average travel time', 'rewards','delay'])
                return f"Your final answer should include this sentence without changing anything: The simulation results are saved at:`{saved_image}` and the log file is saved at:{log_file_path}."
                
            elif algorithm.lower() in support_algorithms_2:
                timestamp = datetime.datetime.now().strftime("%Y_%m_%d-%H_%M_%S")
                
                subprocess.run(["python3", self.base_dir+"run.py", "-a", algorithm,"-tt", timestamp], capture_output=True, text=True, cwd=libsignal_root)

                sub_half = simulator + "_" + algorithm + "/cityflow1x1/0/logger/example_BRF.log"
                path = self.temp_directfeedback + sub_half

                with open(file=path, mode="r") as reader:
                    simulation_result = reader.read()                 
                return f"Your final answer should include this sentence without changing anything: The simulation results are: {simulation_result}."
        
class filterNetwork:
    def __init__(self, figfolder: str) -> None:
        # base_network: 
        self.figfolder = figfolder
        self.store_base = filter_store_base
    
    # The input most likely contains the information related to keywords like: filter network, walk/walkable, bike/bikeable/, railway/railway routes.
    @func_prompt(name='networkFilter',
             description="""
            This tool is used to filter the network by the required categories.
            This tool will also return the file path of a filtered road network with emphasized lanes of interest to provide the final answer. 
            Consider Using this tool if the question is about filtering or obtaining a walkable/bikeable/railway/ network from the base network.
            The output will tell you whether you have finished this command successfully.
             """)

    def embody(self, target: str) -> str:
        # get the target information:
        print("target:")
        print(target)
        # target = target.replace("[", "").replace("]", "").replace(" ", "")
        

        target_path, keyword = target.split(",")
        time.sleep(0.5)
        try:
            net = og.getNetFromFile(target_path, network_types=keyword)
        except Exception as e :
            print("There are some error when filtering")

        time_now = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

        store_info = self.store_base + time_now + "-" + keyword +".png"

        og.osmnet.visualization.saveFig(network=net, picpath=store_info)
        if not os.path.exists(store_info):
            return f"You cannot successfully filter the network."
        
        return f"You have successfully filter the network by type: {target} on the target network. And your final answer should include this sentence without changing anything except for translation: The location of interested {target} is kept at: `{store_info}`."


class constructYardNetwork:
    def __init__(self, figfolder: str) -> None:
        # base_network: 
        self.figfolder = figfolder
        self.store_base = filter_store_base
    
    # The input most likely contains the information related to keywords like: construct yard network, railway/railway network, construct railway yard.
    @func_prompt(name='constructYardNetwork',
            description="""
            This tool is used to construct railway yard network.
            This tool will also return the file path of the constructed yard network.
            Consider using this tool if the question is about filtering or obtaining a railway network from the base network.
            The output will tell you whether you have finished this command successfully.
             """)
    def embody(self, target: str) -> str:
        # get the target information:
        print("target:")
        print(target)
        
        target_path, keyword = target.split(",")
        print("keyword ", keyword)
        # target_path = "/home/local/ASURITE/longchao/Desktop/project/LLM4Traffic/OpenTI/pivotalAgent/Data/download/OSM/IHB_yard.osm"
        time.sleep(0.5)
        try:
            net=orl.get_network_from_file(target_path,POIs=True,check_boundary=True)

        except Exception as e :
            print("There are some error when filtering")

        time_now = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

        store_info = self.store_base + time_now + "-" + keyword +".png"

        orl.show_network(net, store_info)
        orl.save_network(net,output_folder='./csvfile')
        if not os.path.exists(store_info):
            return f"You cannot successfully filter the network."
        
        return f"You have successfully constructed the yard network by type: {target} on the target network. And your final answer should include this sentence without changing anything except for translation: The location of interested {target} is kept at: `{store_info}`."

class generatePolygon:
    # The input most likely contains the information related to keywords like: construct yard network, railway/railway network, construct railway yard.
    @func_prompt(name='generatePolygon',
            description="""
            This tool is used to generate polygon shape for the railway yard network.
            This tool will also return the file path of the constructed polygon file for the yard.
            Consider using this tool if the question is about creating/finding/generating shape of polygon for the railway yard.
            The output will tell you whether you have finished this command successfully.
             """)
    def embody(self, target: str) -> str:

        input_file = target
        # Load the GeoJSON file using json
        with open(input_file, 'r') as f:
            geojson_data = json.load(f)

        # Function to calculate the distance between two points
        def calculate_distance(coord1, coord2):
            return math.sqrt((coord1[0] - coord2[0]) ** 2 + (coord1[1] - coord2[1]) ** 2)

        # Function to find the longest pair of nodes in the polygon (long axis)
        def find_longest_pair(polygon):
            longest_distance = 0
            long_axis_nodes = None
            coords = list(polygon.exterior.coords)
            for i in range(len(coords)):
                for j in range(i + 1, len(coords)):
                    distance = calculate_distance(coords[i], coords[j])
                    if distance > longest_distance:
                        longest_distance = distance
                        long_axis_nodes = (coords[i], coords[j])
            return long_axis_nodes

        # Function to create a minimized polygon with margin 'd' around the MultiLineString or LineString geometry
        def create_envelope_polygon(coordinates, d):
            if not coordinates or not isinstance(coordinates, list):
                return None
            try:
                # Create a LineString or MultiLineString geometry
                lines = MultiLineString(coordinates) if len(coordinates) > 1 else LineString(coordinates[0])

                # Create a convex hull polygon that envelopes the lines
                envelope_polygon = lines.convex_hull

                # Create a buffer around the polygon (to simulate margin 'd')
                buffered_polygon = envelope_polygon.buffer(d / 111320)  # Convert meters to degrees (approx. conversion)

                # Simplify the polygon to reduce the number of points (up to 6 points)
                minimized_polygon = buffered_polygon.simplify(0.01, preserve_topology=True)

                # Ensure the polygon has no more than 6 points
                if len(minimized_polygon.exterior.coords) > 6:
                    minimized_polygon = Polygon(minimized_polygon.exterior.coords[:6])

                return minimized_polygon
            except Exception as e:
                print(f"Error creating polygon: {e}")
                return None


        # Set margin 'd' in meters (e.g., 50 meters)
        d = 50

        # Now, process the geometries to generate the polygons for each yard
        yard_data_with_polygons = []
        for feature in geojson_data['features']:
            try:
                # Extract YARDNAME from properties
                properties = feature.get('properties', {})
                yard_name = properties.get('YARDNAME', 'Unknown')

                # Extract coordinates from geometry (which are MultiLineString or LineString)
                geometry = feature.get('geometry', {})
                coordinates = []
                if geometry.get('type') == 'LineString':
                    coordinates = [geometry.get('coordinates', [])]  # Wrap in list for uniformity
                elif geometry.get('type') == 'MultiLineString':
                    coordinates = geometry.get('coordinates', [])

                # Skip rows with missing or malformed data
                if yard_name == 'Unknown' or not coordinates:
                    continue

                # Generate a minimized polygon if valid coordinates exist
                polygon = create_envelope_polygon(coordinates, d)
                long_axis_nodes=find_longest_pair(polygon)
                if polygon:
                    yard_data_with_polygons.append({
                        'YardName': yard_name,
                        'Polygon_WKT': polygon.wkt,
                        'Main_Direction':long_axis_nodes
                    })
                else:
                    print(f"Invalid geometry for Yard: {yard_name}")
            except Exception as e:
                print(f"Error processing feature: {e}")

        # Convert the data into a DataFrame for easier manipulation
        yard_polygons_df = pd.DataFrame(yard_data_with_polygons)
        yard_polygons_df.index += 1  # Auto-index starting from 1

        # Insert an auto-increment index as the first column
        yard_polygons_df.reset_index(drop=True, inplace=True)
        yard_polygons_df.index += 1  # Auto-index starting from 1
        yard_polygons_df.insert(0, 'ID', yard_polygons_df.index)

        # Remove any empty or invalid columns (if any)
        yard_polygons_df = yard_polygons_df[['ID', 'YardName', 'Polygon_WKT','Main_Direction']]

        # Save the updated DataFrame with the auto-index to a new CSV file
        output_file = os.path.join(os.path.dirname(input_file), "rail_yards_polygons.csv")
        # csv_file_with_autoindex = "input/rail_yards_polygons_IHB_Yard.csv"
        yard_polygons_df.to_csv(output_file, index=False)

class generateBNetwork:
    @func_prompt(name='generateBNetwork',
                 description="""
                 This tool is used to generate a B network from node, link, and polygon CSV files.
                 It returns the file path of the constructed B network.
                 Use this tool to construct a railway network based on provided CSV inputs.
                 """)
    def embody(self, target: str) -> str:
        print("target ", target)
        root_path, keyword = target.split(",")
        node_file = root_path+"/node.csv"
        link_file = root_path+"/link.csv"
        polygon_file = root_path + "/rail_yards_polygons.csv"
        

        # 1. Read the origin network
        ## /home/local/ASURITE/longchao/Desktop/project/LLM4Traffic/OpenTI/RailYardNetwork/
        all_node = pd.read_csv(node_file, low_memory=False)
        all_link = pd.read_csv(link_file, low_memory=False)
        all_yard = pd.read_csv(polygon_file, low_memory=False)
        # all_yard = pd.read_csv("./rail_yards_polygons.csv", low_memory=False)

        all_gate = pd.DataFrame(columns=['gate_node_id', 'belong_yard', 'x_coord', 'y_coord', 
                                        'direction','yard_node_id', 'closest_in', 'closest_out'])
        all_bridge = pd.DataFrame()
        all_yard_node = pd.DataFrame()
        all_yard_link = pd.DataFrame()

        # for key, yard in all_yard.iloc[:10].iterrows():
        for key, yard in all_yard.iterrows():
            try:
                yard_id = yard['ID']
                print(yard_id)
                polygon = wkt.loads(yard['Polygon_WKT'])

                max_x, min_x = max([point[0] for point in polygon.exterior.coords]), min([point[0] for point in polygon.exterior.coords])
                max_y, min_y = max([point[1] for point in polygon.exterior.coords]), min([point[1] for point in polygon.exterior.coords])

                # 2. NODE
                # identify nodes in the rectangle yard zone (roughly)
                yard_node_temp = all_node[all_node['x_coord'] <= max_x]
                yard_node_temp = yard_node_temp[all_node['x_coord'] >= min_x]
                yard_node_temp = yard_node_temp[all_node['y_coord'] <= max_y]
                yard_node_temp = yard_node_temp[all_node['y_coord'] >= min_y]

                # identify nodes in the yard zone (accurately)
                yard_node = pd.DataFrame()
                for key_p, node in yard_node_temp.iterrows():
                    point = Point(node['x_coord'], node['y_coord'])
                    is_within = polygon.contains(point)
                    if (is_within == True):
                        yard_node = pd.concat([yard_node, node.to_frame().T], ignore_index=True)
                yard_node_id_list = yard_node['node_id'].to_list()
                
                # delete nodes from the origin network and generate nodes for A network
                all_node = all_node.merge(yard_node, how='left', indicator=True)
                all_node = all_node[all_node['_merge'] == 'left_only'].drop('_merge', axis=1)

                # create nodes for B network in the yards
                yard_node['yard_id'] = yard_id
                all_yard_node = pd.concat([all_yard_node, yard_node], ignore_index=True)

                # 3. LINK
                # identify links within the yards totally
                yard_link = all_link[all_link['from_node_id'].isin(yard_node_id_list)]
                yard_link = yard_link[all_link['to_node_id'].isin(yard_node_id_list)]

                # delete links from the origin network and generate links for A network
                all_link = all_link.merge(yard_link, how='left', indicator=True)
                all_link = all_link[all_link['_merge'] == 'left_only'].drop('_merge', axis=1)

                # create links for B network in the yards
                yard_link['yard_id'] = yard_id
                all_yard_link = pd.concat([all_yard_link, yard_link], ignore_index=True)

                # 4. BRIDGE
                # identify incoming bridge links and incoming bridge links from the origin network
                in_bridge_link = all_link[all_link['to_node_id'].isin(yard_node_id_list)]
                all_link = all_link.merge(in_bridge_link, how='left', indicator=True)
                all_link = all_link[all_link['_merge'] == 'left_only'].drop('_merge', axis=1)

                # extract incoming bridge links
                in_bridge_link['yard_id'] = yard_id
                in_bridge_link['direction'] = 'in'
                all_bridge = pd.concat([all_bridge, in_bridge_link], ignore_index=True)

                # extract gate nodes with respect to incoming bridge links
                for key, link in in_bridge_link.iterrows():
                    gate_node_id = link['from_node_id']
                    gate_node = all_node[all_node['node_id'] == gate_node_id].squeeze()
                    x_coord = gate_node['x_coord']
                    y_coord = gate_node['y_coord']
                    all_gate.loc[len(all_gate)] = [gate_node['node_id'], yard_id, x_coord, y_coord, 'in', '', '', '']

                # identify outgoing bridge links and incoming bridge links from the origin network
                out_bridge_link = all_link[all_link['from_node_id'].isin(yard_node_id_list)]
                all_link = all_link.merge(out_bridge_link, how='left', indicator=True)
                all_link = all_link[all_link['_merge'] == 'left_only'].drop('_merge', axis=1)

                # extract outgoing bridge links
                out_bridge_link['yard_id'] = yard_id
                out_bridge_link['direction'] = 'out'
                all_bridge = pd.concat([all_bridge, out_bridge_link], ignore_index=True)

                # extract gate nodes with respect to outgoing bridge links
                for key, link in out_bridge_link.iterrows():
                    gate_node_id = link['to_node_id']
                    gate_node = all_node[all_node['node_id'] == gate_node_id].squeeze()
                    x_coord = gate_node['x_coord']
                    y_coord = gate_node['y_coord']
                    all_gate.loc[len(all_gate)] = [gate_node['node_id'], yard_id, x_coord, y_coord, 'out', '', '', '']
                
                # identify gate nodes on the main direction
                coordinates = ast.literal_eval(yard['Main_Direction'])
                start_point, end_point = coordinates[0], coordinates[1] 
                main_direction = LineString([start_point, end_point])
                min_distance_start, min_distance_end = float('inf'), float('inf')
                closest_point_start, closest_point_end = None, None
                for key, gate in all_gate.iterrows():
                    point_gate = Point(gate['x_coord'], gate['y_coord'])

                    dist_start = point_gate.distance(Point(start_point))
                    dist_end = point_gate.distance(Point(end_point))
                    
                    if dist_start < min_distance_start:
                        min_distance_start = dist_start
                        closest_point_start = key
                    if dist_end < min_distance_end:
                        min_distance_end = dist_end
                        closest_point_end = key

                all_gate.at[closest_point_start, 'closest_in'] = 1
                all_gate.at[closest_point_end, 'closest_out'] = 1

            except Exception as e:
                print("ERROR with yard " + str(yard_id))
                continue

        # Define your output folder
        output_folder = './csvfile'

        # Ensure the output folder exists
        os.makedirs(output_folder, exist_ok=True)

        # 5. Merge network
        all_node['resolution_type'] = 'Macro'
        all_yard_node['resolution_type'] = 'Micro'
        all_gate['resolution_type'] = 'Gate'
        cross_resolution_node = pd.concat([all_node, all_yard_node, all_gate], axis = 0) 
        cross_resolution_link = pd.concat([all_link, all_yard_link], axis = 0) 
        cross_resolution_node.to_csv('C_Network_node_testing.csv', index = False)
        cross_resolution_link.to_csv('C_Network_link_testing.csv', index = False)


        # 6. output individual files
        # Save the DataFrames to the specified folder
        all_node.to_csv(os.path.join(output_folder, 'B_Network_node.csv'), index=False)  # B network node
        all_link.to_csv(os.path.join(output_folder, 'B_Network_link.csv'), index=False)  # B network link
        all_gate.to_csv(os.path.join(output_folder, 'gate.csv'), index=False)  # gate node
        all_bridge.to_csv(os.path.join(output_folder, 'bridge.csv'), index=False)  # bridge link
        all_yard_node.to_csv(os.path.join(output_folder, 'B_Network_yard_node.csv'), index=False)  # B network yard node
        all_yard_link.to_csv(os.path.join(output_folder, 'B_Network_yard_link.csv'), index=False)  # B network yard link

        print('finished!')




class generateANetwork:
    @func_prompt(name='generateANetwork',
                 description="""
                 This tool is used to generate a A network from gate and polygon CSV files.
                 It returns the file path of the constructed A network.
                 Use this tool to construct a railway network based on provided CSV inputs.
                 """)
    def embody(self, target: str) -> str:
        # Function to generate unique IDs for nodes and links
        def generate_id():
            return str(uuid.uuid4())


        # Function to calculate the distance between two points
        def calculate_distance(coord1, coord2):
            return math.sqrt((coord1[0] - coord2[0]) ** 2 + (coord1[1] - coord2[1]) ** 2)


        # Function to find the longest pair of nodes in the polygon (long axis)
        def find_longest_pair(polygon):
            longest_distance = 0
            long_axis_nodes = None
            coords = list(polygon.exterior.coords)
            for i in range(len(coords)):
                for j in range(i + 1, len(coords)):
                    distance = calculate_distance(coords[i], coords[j])
                    if distance > longest_distance:
                        longest_distance = distance
                        long_axis_nodes = (coords[i], coords[j])
            return long_axis_nodes


        # Function to assign long_axis_in and long_axis_out based on proximity to gates
        def assign_in_out_nodes(long_axis_nodes, gate_in_coords, gate_out_coords):
            in_node, out_node = long_axis_nodes

            # Find the closest node to gate in
            if gate_in_coords:
                closest_gate_in = min(gate_in_coords, key=lambda gate: calculate_distance(in_node, gate))
                dist_in_node = calculate_distance(in_node, closest_gate_in)
                dist_out_node = calculate_distance(out_node, closest_gate_in)

                if dist_out_node < dist_in_node:
                    in_node, out_node = out_node, in_node  # Swap nodes

            # Ensure that the other node is closest to the gate out
            if gate_out_coords:
                closest_gate_out = min(gate_out_coords, key=lambda gate: calculate_distance(out_node, gate))
                dist_in_node = calculate_distance(in_node, closest_gate_out)
                dist_out_node = calculate_distance(out_node, closest_gate_out)

                if dist_in_node < dist_out_node:
                    in_node, out_node = out_node, in_node  # Swap back if needed

            return in_node, out_node


        # Function to convert meters to degrees (for offsetting siding tracks)
        def meters_to_degrees(offset_meters, y_coord):
            lat_offset = offset_meters / 111320  # 1 degree latitude ≈ 111,320 meters
            lon_offset = offset_meters / (111320 * math.cos(math.radians(y_coord)))  # Adjust longitude based on latitude
            return lon_offset, lat_offset


        # Read rail yard polygon and output gate CSV files
        print("target ", target)
        root_path, keyword = target.split(',')
        gate_file = root_path+"/gate.csv"
        polygon_file = root_path + "/rail_yards_polygons.csv"
        yards_df = pd.read_csv(polygon_file)
        gates_df = pd.read_csv(gate_file)

        # Initialize lists to store nodes and links
        nodes_data = []
        links_data = []

        # Parameter to control the number of siding tracks on each side of the main track
        num_siding_tracks_per_side = 2  # Number of tracks on each side
        track_offset_meters = 150  # Distance between parallel siding tracks in meters
        length_decrease_factor = 0.9  # Factor to decrease the length of each siding track

        # A dictionary to map node IDs to their coordinates for constructing the LINESTRING geometry later
        node_coordinates = {}

        # Process each row in rail_yards_polygons.csv (each yard)
        for index, row in yards_df.iterrows():
            yard_name = row['YardName']
            yard_id = row['ID']
            polygon_wkt = row['Polygon_WKT']

            # Load the polygon geometry from WKT
            polygon = load_wkt(polygon_wkt)

            # Find the longest pair of nodes (long axis)
            long_axis_nodes = find_longest_pair(polygon)
            main_track_length = calculate_distance(long_axis_nodes[0], long_axis_nodes[1])

            # Find gates in output_gate.csv for the current yard
            yard_gates = gates_df[gates_df['belong_yard'] == yard_id]
            gate_in_coords = [(row['x_coord'], row['y_coord']) for index, row in yard_gates[yard_gates['direction'] == 'in'].iterrows()]
            gate_out_coords = [(row['x_coord'], row['y_coord']) for index, row in yard_gates[yard_gates['direction'] == 'out'].iterrows()]

            # Assign long_axis_in_node and long_axis_out_node based on proximity to gates
            long_axis_in_node, long_axis_out_node = assign_in_out_nodes(long_axis_nodes, gate_in_coords, gate_out_coords)

            # Add the long_axis_in and long_axis_out nodes to the node.csv
            long_axis_in_node_id = generate_id()
            long_axis_out_node_id = generate_id()

            nodes_data.append({
                'node_id': long_axis_in_node_id,
                'yard_name': yard_name,
                'node_type': 'long_axis_in',
                'x_coord': long_axis_in_node[0],
                'y_coord': long_axis_in_node[1]
            })

            nodes_data.append({
                'node_id': long_axis_out_node_id,
                'yard_name': yard_name,
                'node_type': 'long_axis_out',
                'x_coord': long_axis_out_node[0],
                'y_coord': long_axis_out_node[1]
            })

            # Store coordinates for constructing LINESTRING geometries later
            node_coordinates[long_axis_in_node_id] = long_axis_in_node
            node_coordinates[long_axis_out_node_id] = long_axis_out_node

            # Add the long axis link to link.csv
            links_data.append({
                'from_node_id': long_axis_in_node_id,
                'to_node_id': long_axis_out_node_id,
                'link_type': 'main_track',
                'geometry': f"LINESTRING({long_axis_in_node[0]} {long_axis_in_node[1]}, {long_axis_out_node[0]} {long_axis_out_node[1]})"
            })

            # Add gate_in and gate_out points to node.csv
            for gate_in in gate_in_coords:
                gate_in_node_id = generate_id()
                nodes_data.append({
                    'node_id': gate_in_node_id,
                    'yard_name': yard_name,
                    'node_type': 'gate_in',
                    'x_coord': gate_in[0],
                    'y_coord': gate_in[1]
                })
                node_coordinates[gate_in_node_id] = gate_in  # Store coordinates

                links_data.append({
                    'from_node_id': gate_in_node_id,
                    'to_node_id': long_axis_in_node_id,
                    'link_type': 'bridge_in',
                    'geometry': f"LINESTRING({gate_in[0]} {gate_in[1]}, {long_axis_in_node[0]} {long_axis_in_node[1]})"
                })

            for gate_out in gate_out_coords:
                gate_out_node_id = generate_id()
                nodes_data.append({
                    'node_id': gate_out_node_id,
                    'yard_name': yard_name,
                    'node_type': 'gate_out',
                    'x_coord': gate_out[0],
                    'y_coord': gate_out[1]
                })
                node_coordinates[gate_out_node_id] = gate_out  # Store coordinates

                links_data.append({
                    'from_node_id': gate_out_node_id,
                    'to_node_id': long_axis_out_node_id,
                    'link_type': 'bridge_out',
                    'geometry': f"LINESTRING({gate_out[0]} {gate_out[1]}, {long_axis_out_node[0]} {long_axis_out_node[1]})"
                })

            # Generate siding tracks (parallel to the main track) on both sides
            lon_offset, lat_offset = meters_to_degrees(track_offset_meters, long_axis_in_node[1])

            for i in range(1, num_siding_tracks_per_side + 1):
                # Decrease length of the siding track outward
                track_length = main_track_length * (length_decrease_factor ** i)
                length_ratio = track_length / main_track_length

                # Compute new coordinates based on length_ratio, keeping it parallel to the main track
                siding_in_node_pos = (long_axis_in_node[0], long_axis_in_node[1] + i * lat_offset)
                siding_out_node_pos = (long_axis_out_node[0], long_axis_out_node[1] + i * lat_offset)

                siding_in_node_neg = (long_axis_in_node[0], long_axis_in_node[1] - i * lat_offset)
                siding_out_node_neg = (long_axis_out_node[0], long_axis_out_node[1] - i * lat_offset)

                siding_in_node_pos_id = generate_id()
                siding_in_node_neg_id = generate_id()
                siding_out_node_pos_id = generate_id()
                siding_out_node_neg_id = generate_id()

                # Add nodes and links for positive and negative siding tracks
                nodes_data.append({
                    'node_id': siding_in_node_pos_id,
                    'yard_name': yard_name,
                    'node_type': 'track',
                    'x_coord': siding_in_node_pos[0],
                    'y_coord': siding_in_node_pos[1]
                })

                nodes_data.append({
                    'node_id': siding_out_node_pos_id,
                    'yard_name': yard_name,
                    'node_type': 'track',
                    'x_coord': siding_out_node_pos[0],
                    'y_coord': siding_out_node_pos[1]
                })

                nodes_data.append({
                    'node_id': siding_in_node_neg_id,
                    'yard_name': yard_name,
                    'node_type': 'track',
                    'x_coord': siding_in_node_neg[0],
                    'y_coord': siding_in_node_neg[1]
                })

                nodes_data.append({
                    'node_id': siding_out_node_neg_id,
                    'yard_name': yard_name,
                    'node_type': 'track',
                    'x_coord': siding_out_node_neg[0],
                    'y_coord': siding_out_node_neg[1]
                })

                # Store coordinates for LINESTRING construction
                node_coordinates[siding_in_node_pos_id] = siding_in_node_pos
                node_coordinates[siding_out_node_pos_id] = siding_out_node_pos
                node_coordinates[siding_in_node_neg_id] = siding_in_node_neg
                node_coordinates[siding_out_node_neg_id] = siding_out_node_neg

                # Add positive siding track links to link.csv (parallel to the main track)
                links_data.append({
                    'from_node_id': siding_in_node_pos_id,
                    'to_node_id': siding_out_node_pos_id,
                    'link_type': 'siding_track',
                    'geometry': f"LINESTRING({siding_in_node_pos[0]} {siding_in_node_pos[1]}, {siding_out_node_pos[0]} {siding_out_node_pos[1]})"
                })

                # Add negative siding track links to link.csv (parallel to the main track)
                links_data.append({
                    'from_node_id': siding_in_node_neg_id,
                    'to_node_id': siding_out_node_neg_id,
                    'link_type': 'siding_track',
                    'geometry': f"LINESTRING({siding_in_node_neg[0]} {siding_in_node_neg[1]}, {siding_out_node_neg[0]} {siding_out_node_neg[1]})"
                })

                # Connect siding track in_nodes to the long_axis_in_node
                links_data.append({
                    'from_node_id': siding_in_node_pos_id,
                    'to_node_id': long_axis_in_node_id,
                    'link_type': 'track',
                    'geometry': f"LINESTRING({siding_in_node_pos[0]} {siding_in_node_pos[1]}, {long_axis_in_node[0]} {long_axis_in_node[1]})"
                })

                links_data.append({
                    'from_node_id': siding_in_node_neg_id,
                    'to_node_id': long_axis_in_node_id,
                    'link_type': 'track',
                    'geometry': f"LINESTRING({siding_in_node_neg[0]} {siding_in_node_neg[1]}, {long_axis_in_node[0]} {long_axis_in_node[1]})"
                })

                # Connect siding track out_nodes to the long_axis_out_node
                links_data.append({
                    'from_node_id': siding_out_node_pos_id,
                    'to_node_id': long_axis_out_node_id,
                    'link_type': 'track',
                    'geometry': f"LINESTRING({siding_out_node_pos[0]} {siding_out_node_pos[1]}, {long_axis_out_node[0]} {long_axis_out_node[1]})"
                })

                links_data.append({
                    'from_node_id': siding_out_node_neg_id,
                    'to_node_id': long_axis_out_node_id,
                    'link_type': 'track',
                    'geometry': f"LINESTRING({siding_out_node_neg[0]} {siding_out_node_neg[1]}, {long_axis_out_node[0]} {long_axis_out_node[1]})"
                })

        # Convert node and link data to DataFrames and save to CSV
        nodes_df = pd.DataFrame(nodes_data)
        links_df = pd.DataFrame(links_data)

        # Define your output folder
        output_folder = './csvfile'

        # Ensure the output folder exists
        os.makedirs(output_folder, exist_ok=True)

        # Save to CSV, now with the geometry column in the links.csv file
        nodes_df.to_csv(os.path.join(output_folder,'A_Network_node.csv'), index=False)
        links_df.to_csv(os.path.join(output_folder,'A_Network_link.csv'), index=False)

        print("Nodes and links with geometry have been saved to A_Network_node.csv and A_Network_link.csv")


class generateDemand:

    def __init__(self, save: str) -> None:
        self.demand_save = save
    
    @func_prompt(name="generateDemand", 
             description="""This tool is used for generating demand based on the osm data, it leverages the package named "grid2demand".
             Consider using this tool when asked to generate traffic demand based on the downloaded map/osm file. And you will need the a string of path to the .osm file.
             The output will tell you whether you have finished this command successfully.
             """)
    
    def embody(self, target: str) -> str:
        import grid2demand as gd

        print("demand target:{}".format(target))
        name = target.replace(".osm", "").split("/")[-1]
        dump_dir = "./Data/demand/"+name
        net = og.getNetFromFile(target, network_types=('walk', 'auto') , POI=True, default_lanes=True, default_speed=True)    
        og.connectPOIWithNet(net)
        og.generateNodeActivityInfo(net)
        if not os.path.exists(dump_dir):
            os.makedirs(dump_dir)
        og.outputNetToCSV(net, output_folder=dump_dir)
        gd.read_network(dump_dir)

        # Initialize a GRID2DEMAND object
        gd = GRID2DEMAND(dump_dir)

        # Load node and poi data from input directory
        node_dict, poi_dict = gd.load_network.values()

        # Generate zone dictionary from node dictionary by specifying number of x blocks and y blocks
        zone_dict = gd.net2zone(node_dict, num_x_blocks=20, num_y_blocks=20)

        # synchronize geometry info between zone, node and poi
        # add zone_id to node and poi dictionaries
        # also add node_list and poi_list to zone dictionary
        updated_dict = gd.sync_geometry_between_zone_and_node_poi(zone_dict, node_dict, poi_dict)
        zone_dict_update, node_dict_update, poi_dict_update = updated_dict.values()

        # Calculate zone-to-zone od distance matrix
        zone_od_distance_matrix = gd.calc_zone_od_distance_matrix(zone_dict_update)

        # Generate poi trip rate for each poi
        poi_trip_rate = gd.gen_poi_trip_rate(poi_dict_update)

        # Generate node production attraction for each node based on poi_trip_rate
        node_prod_attr = gd.gen_node_prod_attr(node_dict_update, poi_trip_rate)

        #Calculate zone production and attraction based on node production and attraction
        zone_prod_attr = gd.calc_zone_prod_attr(node_prod_attr, zone_dict_update)

        #Run gravity model to generate agent-based demand
        df_demand = gd.run_gravity_model(zone_prod_attr, zone_od_distance_matrix)

        # Generate agent-based demand
        df_agent = gd.gen_agent_based_demand(node_prod_attr, zone_prod_attr, df_demand=df_demand)

        # You can also view and edit the package setting by using gd.pkg_settings
        print(gd.pkg_settings)

        #Output demand, agent, zone, zone_od_dist_table, zone_od_dist_matrix files to output directory
        gd.save_demand
        df_demand.to_csv(dump_dir+"/demand2.csv")
        gd.save_agent
        gd.save_zone
        gd.save_zone_od_dist_table
        gd.save_zone_od_dist_matrix

        try:
            subprocess.run(['xdg-open', dump_dir], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Failed to open the file: {e}")


        return f"You have successfully generated the demand files. And your final answer should include this sentence without changing anything except for translation: The generated demand is kept at: {dump_dir}."


class simulateOnDLSim:
    def __init__(self, demand_path: str) -> None:
        self.demand_path = demand_path
        self.simulate_path = "./AugmentUtils/simulation/"

    @func_prompt(name="simulateOnDLSim", description="""
    This tool is used for simulating on the DLSim multi-resolution traffic simulator. 
    Please consider using this tool when asked to run simulation on DLSim simulator given a demand path. Please try to provide log info in bullet list or raw data, even though it is a large file.
    """)

    def embody(self, target: str) -> str:
        target = "./AugmentUtils/simulation/simulate/"
        print("DLSim target:{}".format(target))
        try: 
            subprocess.run(["python3", self.simulate_path+"simulate.py"])
        except Exception as e:
            print(e)
        time_now = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

        with open(config_path, 'r') as file:
            base_path = yaml.safe_load(file)["DLSim"]["simulate_base"]
            
        save_to = base_path + time_now + "-log.txt"

        log_file = base_path + "log.txt"
        with open(file=log_file, mode="r") as content:
            log_data = content.read()
            keep = log_data.split("Step 1.7:")[-1]
            output = "Step 1.7:" + keep
        
        with open(file=save_to, mode="w+") as new_file: 
            new_file.write(log_data)

        print("output:")
        print(output)

        try:
            subprocess.run(['xdg-open', save_to], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Failed to open the file: {e}")

        return f"You have successfully simulated on DLSim. And your final answer should include this sentence without changing anything: The simulation process and logs are saved below: `{save_to}`, the key steps as: {output}."
        
    

class visualizeDemand:
    def __init__(self, demand_path: str) -> None:
        self.demand_path = demand_path

    @func_prompt(name="visualizeDemand", description="""
    This tool is used for visualizing the demand file generated automatically. 
    Please consider using this tool when asked to present a visualization on the generated demand file.
    The output will tell you whether you have finished this command successfully.
    """)

    def embody(self, target: str) -> str:

        print("DLSim target:{}".format(target))

        with open(config_path, 'r') as file:
            example_root = yaml.safe_load(file)["Demand"]["example"]
        
        base_path = example_root

        return f"You have successfully visualized the traffic demand information. And your final answer should include this sentence without changing anything: The traffic demand information at ASUtempe as below: `{base_path}`."
    

class log_analyzer:
    def __init__(self) -> None:
        pass
    
    @func_prompt(name="logAnalyzer", description="""
    This is the tool used to analyze the log files and provide comparison, if you are asked to analyze any given loaction files, please try your best to find the data, and provide logical and rational understanding on it.
    The output will tell you whether you have finished this command successfully.
    """)

    def embody(self, target: str):
        try:
            #process the input target information and provide the acceptable form
            path = target 
            with open(file=path, mode="r") as content:
                log_data = content.read()
            return f"You have the content read now, please provide your own understanding and provide the explained result of it, the content is: "+ log_data
        except Exception as e:
            return "Your final answer should include this sentence without changing anything: The path you provided is not valid, please examine again or it is also acceptable to pass a absolute path, thank you!"
        

class oneClickSUMO:
    def __init__(self) -> None:
        self.pipline_path = "/pivotalAgent/AugmentUtils/one_click_SUMO/"
    @func_prompt(name="oneClickSUMOSimulation", description="""
    This is the tool used to execute one click run SUMO simulation, if you are asked to run sumo simulation with one click, please consider using this tool. Before using, please check for an osm file path.
    The output will tell you whether you have finished this command successfully.
    """)

    def embody(self, target: str):
        # try:
            osm_file_path = target
            import os
            current_path = os.getcwd()
            print(os.getcwd())
            print(self.pipline_path)

            if osm_file_path == None: return "Your final answer should include this sentence without changing anything: You are expected to input your osm file path (absolute path) to run this simulation!"
            # subprocess.run(["python3", current_path+ self.pipline_path + "simulate_pipline.py", "-f", target], capture_output=True, text=True, cwd=current_path+ self.pipline_path)
            print("target:")
            print(target)
            target = target.replace("`", "")
            print(target)
            # file_name = os.path.basename(target)
            with open(config_path, 'r') as file:
                sumo_pipline = yaml.safe_load(file)["SUMO_PIPLINE"]
            subprocess.run(["python3", sumo_pipline, "-f", target], capture_output=True, text=True, cwd=current_path+ self.pipline_path)

            return "Successfully executed SUMO simulation with random sampled demand information."
        
        # except Exception as e:
            # return "Your final answer should include this sentence without changing anything: The execution happened with error, please examine the path or format, thank you!"
        


class odMatrixTest:
    def __init__(self) -> None:
        pass

    @func_prompt(name="ODMatrixTest", description="""
    This is the tool used to conduct an OD test (Origin-Destination Test). If user asked about running the OD test, or origin-destination matrix optimization, please consider this approach.
    The output will tell you whether you have finished this command successfully.
    """)

    def embody(self, target: str):
        try:
            import os
            current_path = os.getcwd() + "/pivotalAgent/AugmentUtils/od_optimize/simulation_test/"

            subprocess.run(["python3", current_path + "pickout_test_direct.py"], capture_output=True, text=True, cwd=current_path)
            return "The simulation is running in backend process, please wait (estimate time: 8min...)"
        
        except Exception as e:
            return "Your final answer should include this sentence without changing anything: The execution happened with error, please examine the path or format, thank you!"

import pandas as pd

def vialize_od_results(file_path):

    df = pd.read_csv(file_path, na_values='', keep_default_na=False).fillna(0)

    volume = df['volume'].tolist()
    obs_count = df['obs_count'].tolist()


    len_vol = len(volume)
    len_obs = len(obs_count)

    if len_vol == len_obs:
        for i in range(len_vol):
            flag = 0
            if volume[i] == 0:
                flag = 1
            if obs_count[i] == 0:
                flag = 1
            if flag == 1:
                volume[i] = 0
                obs_count[i] = 0



    # combined = [0 if v == 0 or o == 0 else 1 for v, o in zip(volume, obs_count)]

    #  plot the bar chart
    indexes = range(len(volume))

    # Set up the matplotlib figure and axes
    fig, ax = plt.subplots()

    # Add bars for 'volume' and 'obs_count'. The indexes are shifted for each bar for clarity.
    ax.bar([i - 0.2 for i in indexes], volume, width=0.4, label='Ground Truth Volume', align='center')
    ax.bar([i + 0.2 for i in indexes], obs_count, width=0.4, label='OD Matrix Count', align='center')

    # Add bars for 'combined' list to reflect positions with 0 in either 'volume' or 'obs_count'
    # ax.bar(indexes, combined, width=0.4, label='Combined', color='r', alpha=0.5, align='center')

    # Labeling
    ax.set_xlabel('Share Index')
    ax.set_ylabel('Values')
    ax.set_title('Volume and Obs Count Comparison')
    ax.legend()

    # Show the plot
    plt.show()


class visualizeODResult:
    def __init__(self) -> None:
        pass

    @func_prompt(name="visualizeODResult", description="""
    This is the tool used to visualize the OD matrix estimation results from DTALite simulator, If user ask about visualize the OD results of specific file path of DTALite execution, you should consider using this tool.
    The output will tell you whether you have finished this command successfully.              
    """)

    def embody(self, target: str):
        try:
            if os.path.exists(target):
                vialize_od_results(target)
            else:
                return "There is no such a file after execution of DTALite, please check again."

        except Exception as e:
            return "Your final answer should include this sentence without changing anything: The execution happened with error, please examine the path or format, thank you!"
        