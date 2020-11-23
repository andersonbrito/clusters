# coding=utf-8
import pycountry_convert as pyCountry
import pandas as pd
import argparse
from bs4 import BeautifulSoup as BS
import pycountry
from matplotlib import cm
import numpy as np


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Generate ordered colour file for nextstrain build",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", required=True, help="Reformatted nextstrain metadata file")
    parser.add_argument("--coordinates", required=True,  help="TSV coordinates file being used in the build")
    parser.add_argument("--geoscheme", required=True, help="XML file with geographic scheme")
    parser.add_argument("--grid", required=True, help="HTML file with HEX colour matrices")
    parser.add_argument("--columns", nargs='+', type=str,   help="list of columns with geographic information")
    parser.add_argument("--output", required=True, help="TSV file containing ordered HEX colours based on locations")
    args = parser.parse_args()

    metadata = args.metadata
    geoscheme = args.geoscheme
    coordinates = args.coordinates
    grid = args.grid
    columns = args.columns
    output = args.output


    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov/ncov_nfl/nextstrain/batchXX_newcolours/pre-analyses/'
    # metadata = path + 'metadata_geo.tsv'
    # coordinates = path + 'latlongs.tsv'
    # geoscheme = path + 'geoscheme.tsv'
    # grid = path + 'colour_grid.html'
    # columns = ['region_exposure', 'country_exposure', 'division_exposure', 'location']
    # output = path + 'colors.tsv'


    # pre-determined HEX colours and hues
    force_colour = {}
    force_hue = {'North America': 0}

    # content to be exported as final result
    latlongs = {trait: {} for trait in columns}

    # get ISO alpha3 country codes
    isos = {}
    def get_iso(country):
        global isos
        if country not in isos.keys():
            try:
                isoCode = pyCountry.country_name_to_country_alpha3(country, cn_name_format="default")
                isos[country] = isoCode
            except:
                try:
                    isoCode = pycountry.countries.search_fuzzy(country)[0].alpha_3
                    isos[country] = isoCode
                except:
                    isos[country] = ''
        return isos[country]


    # extract coordinates from latlongs file for sorting places by latitude
    for line in open(coordinates).readlines():
        if not line.startswith('\n'):
            try:
                trait, place, lat, long = line.strip().split('\t')
                if trait in latlongs.keys():
                    entry = {place: (str(lat), str(long))}
                    latlongs[trait].update(entry) # save as pre-existing result
                else:
                    print('### WARNING! ' + trait + ' is not among the pre-selected columns with geographic information!')
            except:
                pass



    ''' REORDER LOCATIONS FOR LEGEND FORMATTING '''
    print('\n### Parsing metadata...\n')
    # open metadata file as dataframe
    dfN = pd.read_csv(metadata, encoding='utf-8', sep='\t', dtype=str)

    # set US-based geoscheme
    scheme_list = open(geoscheme, "r").readlines()[1:]
    us_geoscheme = {}
    for line in scheme_list:
        if not line.startswith('\n'):
            type = line.split('\t')[0]
            if type == 'country_exposure':
                region = line.split('\t')[2]
                states = [state.strip() for state in line.split('\t')[5].split(',')]
                for state in states:
                    if state not in us_geoscheme:
                        us_geoscheme[state] = region

    # collect facility data
    list_facilities = set([facility for facility in dfN['category'].values.tolist() if 'Facility' in facility])
    usregion_facilities = {}
    for fac in list_facilities:
        state = dfN.loc[dfN['category'] == fac, 'division_exposure'].to_list()[0]
        if state in us_geoscheme:
            region = us_geoscheme[state]
            if region not in usregion_facilities:
                usregion_facilities[region] = [fac]
            else:
                usregion_facilities[region] += [fac]

    # keep only geographic data columns
    dfN = dfN[['region_exposure', 'country_exposure', 'division_exposure', 'location']]

    ordered_regions = {}
    dcountries = {}
    places = []
    pinpoints = [dfN[trait].values.tolist() for trait in columns]
    print('\tStep 1. Setting colours for \'' + columns[0] + '\' level...')
    for region_index in latlongs[columns[0]]:
        for address in zip(*pinpoints):
            address = list(address)
            region, country, division, location = [trait for trait in columns]
            region_name = address[0]
            country_name = address[1]
            division_name = address[2]
            location_name = address[3]
            places.append(address)
            # print(address)
            if region_index == region_name:
                if 'subcontinent' not in ordered_regions.keys():
                    ordered_regions['subcontinent'] = {}
                if region_name not in ordered_regions['subcontinent'].keys():
                    ordered_regions['subcontinent'].update({region_name: latlongs[region][region_name]})

                ### STORE COUNTRIES WITH REGIONS AS KEYS
                if country_name in latlongs[country].keys():
                    if region_name not in dcountries.keys():
                        dcountries[region_name] = {country_name: latlongs[country][country_name]}
                    else:
                        if country_name not in dcountries[region_name].keys():
                            dcountries[region_name].update({country_name: latlongs[country][country_name]})

    # sort division entries based on sorted country entries
    print('\tStep 2. Setting colours for \'' + columns[1] + '\' level...')
    ordered_countries = {}
    for region, countries in dcountries.items():
        ordered_countries[region] = {k: v for k, v in sorted(countries.items(), key=lambda item: item[1])}

    print('\tStep 3. Setting colours for \'' + columns[2] + '\' level...')
    ddivisions = {}
    for country_index in [key for dict_ in ordered_countries.values() for key in dict_]:
        for address in places:
            region, country, division, location = [trait for trait in columns]
            country_name = address[1]
            division_name = address[2]
            location_name = address[3]

            if country_index == country_name:
                ### STORE DIVISIONS WITH COUNTRIES AS KEYS
                if division_name in latlongs[division].keys():
                    if country_name not in ddivisions.keys():
                        ddivisions[country_name] = {division_name: latlongs[division][division_name]}
                    else:
                        if division_name not in ddivisions[country_name].keys():
                            ddivisions[country_name].update({division_name: latlongs[division][division_name]})


    # sort division entries based on sorted country entries
    ordered_divisions = {}
    for country, divisions in ddivisions.items():
        ordered_divisions[country] = {k: v for k, v in sorted(divisions.items(), key=lambda item: item[1])}

    dlocations = {}
    for division_index in [key for dict_ in ordered_divisions.values() for key in dict_]:
        for address in places:
            address = list(address)
            region, country, division, location = [trait for trait in columns]
            division_name = address[2]
            location_name = address[3]

            if division_index == division_name:
                ### STORE LOCATIONS WITH DIVISIONS AS KEYS
                if location_name in latlongs[location].keys():
                    if division_name not in dlocations.keys():
                        dlocations[division_name] = {location_name: latlongs[location][location_name]}
                    else:
                        if location_name not in dlocations[division_name].keys():
                            dlocations[division_name].update({location_name: latlongs[location][location_name]})

    print('\tStep 4. Setting colours for \'' + columns[3] + '\' level...')

    # sort locations entries based on sorted division entries
    ordered_locations = {}
    for division, locations in dlocations.items():
        ordered_locations[division] = {k: v for k, v in sorted(locations.items(), key=lambda item: item[1])}



    ''' CONVERT RGB TO HEX COLOURS '''

    # original source: https://bsou.io/posts/color-gradients-with-python

    # convert colour codes
    def RGB_to_hex(RGB):
        ''' [255,255,255] -> "#FFFFFF" '''
        # Components need to be integers for hex to make sense
        RGB = [int(x) for x in RGB]
        return "#" + "".join(["0{0:x}".format(v) if v < 16 else
                              "{0:x}".format(v) for v in RGB])

    def hex_to_RGB(hex):
        ''' "#FFFFFF" -> [255,255,255] '''
        # Pass 16 to the integer function for change of base
        return [int(hex[i:i + 2], 16) for i in range(1, 6, 2)]


    def color_dict(gradient):
        ''' Takes in a list of RGB sub-lists and returns dictionary of
          colors in RGB and hex form for use in a graphing function
          defined later on '''
        return [RGB_to_hex(RGB) for RGB in gradient]

    # create gradient
    def linear_gradient(start_hex, finish_hex, n):
        ''' returns a gradient list of (n) colors between
          two hex colors. start_hex and finish_hex
          should be the full six-digit color string,
          ("#FFFFFF") '''
        # Starting and ending colors in RGB form
        s = hex_to_RGB(start_hex)
        f = hex_to_RGB(finish_hex)
        # Initilize a list of the output colors with the starting color
        RGB_list = [s]
        # Calcuate a color at each evenly spaced value of t from 1 to n
        n = n
        for t in range(1, n):
            # Interpolate RGB vector for color at the current value of t
            curr_vector = [
                int(s[j] + (float(t) / (n - 1)) * (f[j] - s[j]))
                for j in range(3)]
            # Add it to list of output colors
            RGB_list.append(curr_vector)
        RGB_list = RGB_list
        return color_dict(RGB_list)




    ''' IMPORT GEOSCHEME '''

    # xml = BS(open(geoscheme, "r").read(), 'xml')
    # levels = xml.find('levels')

    sampled_region = [key for dict_ in ordered_regions.values() for key in dict_]
    geodata = {}

    for line in scheme_list:
        if not line.startswith('\n'):
            type = line.split('\t')[0]
            if type == 'region_exposure':
                continent = line.split('\t')[1]
                region = line.split('\t')[2]
                if region in sampled_region:
                    if continent not in geodata.keys():
                        geodata[continent] = [region]
                    else:
                        geodata[continent] += [region]


    ''' IMPORT COLOUR SCHEME '''
    print('\tStep 5. Setting special custom colours...')
    print('\n### Generating colour scheme...\n')
    html = BS(open(grid, "r").read(), 'html.parser')
    tables = html.find_all('table')

    # for colour_name in colour_scale.keys():
    limits = {'dark': (60, 40), 'light': (100, 80)}  # define saturation and luminance, max and min, respectively
    hue_to_hex = {}
    for table in tables:
        string_head = str(table.caption)
        if string_head == 'None':
            continue
        else:
            hue_value = int(string_head.split('"')[1])

        lux = table.tbody.find_all('tr')
        lum_value = 90  # brightest colour
        hexdark = ''
        hexligth = ''
        for row in lux:
            sat_value = 10  # unsaturated colour
            for cell in row.find_all('td'):
                # pass
                if sat_value == limits['dark'][0] and lum_value == limits['dark'][1]:
                    hexdark = cell.text.strip()
                if sat_value == limits['light'][0] and lum_value == limits['light'][1]:
                    hexligth = cell.text.strip()
                hex = (hexdark, hexligth)
                sat_value += 10
            lum_value -= 10
        hex = (hexdark, hexligth)
        hue_to_hex[hue_value] = hex


    colour_scale = {'magenta': [320], 'purple': [310, 300, 290, 280, 270, 260],
                    'blue': [250, 240, 230, 220], 'cyan': [210, 200, 190, 180], 'turquoise': [170, 160, 150],
                    'green': [140, 130, 120], 'yellowgreen': [110, 100, 90, 80, 70],
                    'yellow': [60, 50, 40], 'orange': [30, 20], 'red': [10, 0]}

    continent_hues = {'Oceania': colour_scale['magenta'], 'Asia': colour_scale['purple'],
                      'Europe': colour_scale['blue'] + colour_scale['cyan'], 'Africa': colour_scale['yellowgreen'],
                      'America': colour_scale['yellow'] + colour_scale['orange'] + colour_scale['red']}

    colour_wheel = {}
    palette = {}
    for area, subareas in geodata.items():
        num_subareas = len(subareas)
        hues = len(continent_hues[area])
        # print(area, subareas)
        for position, subarea in zip([int(x) for x in np.linspace(0, int(hues), num_subareas, endpoint=False)], subareas):
            if subarea not in palette.keys():
                hue = continent_hues[area][position] # colour picker
                palette[subarea] = hue
                # print(subarea, hue)

    # print(sampled_region)
    # print(palette)
    for region in sampled_region:
        if region in force_hue:
            colour_wheel[region] = hue_to_hex[force_hue[region]]
        else:
            colour_wheel[region] = hue_to_hex[palette[region]]


    results = {trait: {} for trait in columns}


    ''' APPLY SAME HUE FOR MEMBERS OF THE SAME SUB-CONTINENT '''

    # assign countries to regions
    country_colours = {}
    reference_countries = {}
    for region, members in ordered_countries.items():
        countries = list(members.keys())
        for country in countries:
            reference_countries[country] = region
            country_colours[region] = countries


    # assign divisions to countries
    division_colours = {}
    reference_divisions = {}
    for country, members in ordered_divisions.items():
        divisions = list(members.keys())
        for division in divisions:
            reference_divisions[division] = reference_countries[country]
            if reference_countries[country] not in division_colours.keys():
                division_colours[reference_countries[country]] = divisions
            else:
                if division not in division_colours[reference_countries[country]]:
                    division_colours[reference_countries[country]].append(division)


    # assign locations to divisions
    location_colours = {}
    reference_locations = {}
    for division, members in ordered_locations.items():
        locations = list(members.keys())
        for location in locations:
            if reference_divisions[division] not in location_colours.keys():
                location_colours[reference_divisions[division]] = locations
            else:
                if location not in location_colours[reference_divisions[division]]:
                    location_colours[reference_divisions[division]].append(location)

    ''' CREATE COLOUR GRADIENT '''
    # define gradients for regions
    for continent, regions in geodata.items():
        hex_limits = []
        for region in sampled_region:
            if region in regions:
                if continent == 'America':
                    hex_limits += [list(colour_wheel[region])[0]]
                else:
                    hex_limits += list(colour_wheel[region])

        start, end = hex_limits[0], hex_limits[-1]
        if len(regions) == 1:
            gradient = linear_gradient(start, end, 4)
            gradient = [list(gradient)[2]]
        else:
            gradient = linear_gradient(start, end, len(regions))

        for region, colour in zip(regions, gradient):
            print('region_exposure', region, colour)
            results['region_exposure'].update({region: colour})

    # define gradients for country
    for hue, countries in country_colours.items():
        start, end = colour_wheel[hue]
        if len(countries) == 1:
            gradient = linear_gradient(start, end, 4)
            gradient = [list(gradient)[2]]
        else:
            gradient = linear_gradient(start, end, len(countries))
        for country, colour in zip(countries, gradient):
            print('country_exposure', country, colour)
            results['country_exposure'].update({country: colour})

    # define gradients for divisions
    for hue, divisions in division_colours.items():
        start, end = colour_wheel[hue]
        if len(divisions) == 1:
            gradient = linear_gradient(start, end, 3)
            gradient = [list(gradient)[1]]
        else:
            gradient = linear_gradient(start, end, len(divisions))
        for division, colour in zip(divisions, gradient):
            print('division_exposure', division, colour)
            results['division_exposure'].update({division: colour})

    # define gradients for locations
    for hue, locations in location_colours.items():
        start, end = colour_wheel[hue]
        if len(locations) == 1:
            gradient = linear_gradient(start, end, 3)
            gradient = [list(gradient)[1]]
        else:
            gradient = linear_gradient(start, end, len(locations))
        for location, colour in zip(locations, gradient):
            print('location', location, colour)
            results['location'].update({location: colour})


    # special colouring
    categories = {'International': '#BFBFBF', 'USA': '#A6A6A6', 'Community':'#1d3026'}
    results['category'] = {}
    usregion_hues = {
        'USA-Northeast': 220, # blues
        'USA-Midwest': 280, # purples
        'USA-South': 0, # reds
        'USA-West': 40 # yellows
        }

    for us_region, hue in usregion_hues.items():
        start, end = hue_to_hex[hue]
        if us_region in usregion_facilities:
            facilities = sorted(usregion_facilities[us_region])
            gradient = linear_gradient(start, end, len(facilities))
            # print(us_region, hue, facilities, gradient)
            for facility, colour in zip(facilities, gradient):
                categories[facility] = colour

    for cat, hex in categories.items():
        results['category'].update({cat: hex})
        print('category', cat, hex)


    ''' EXPORT COLOUR FILE '''
    with open(output, 'w') as outfile:
        for trait, entries in results.items():
            for item, hexcolour in entries.items():
                if item in force_colour and trait not in ['location']:
                    hexcolour = force_colour[item]
                    print('* ' + item + ' is hardcode with the colour ' + hexcolour)
                line = "{}\t{}\t{}\n".format(trait, item, hexcolour.upper())
                outfile.write(line)
            outfile.write('\n')
    print('\nOrdered colour file successfully created!\n')

