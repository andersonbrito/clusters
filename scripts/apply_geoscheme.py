import pycountry_convert as pyCountry
import pycountry
import pandas as pd
import argparse
from uszipcode import SearchEngine


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Reformat metadata file by adding column with subcontinental regions based on the UN geo-scheme",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--metadata", required=True, help="Nextstrain metadata file")
    parser.add_argument("--geoscheme", required=True, help="XML file with geographic classifications")
    parser.add_argument("--output", required=True, help="Updated metadata file")
    args = parser.parse_args()

    metadata = args.metadata
    geoscheme = args.geoscheme
    output = args.output

    # path = "/Users/anderson/GLab Dropbox/Anderson Brito/projects/ncov_nfl/nextstrain/batch01_20201012e/pre-analyses/"
    # metadata = path + 'metadata_filtered.tsv'
    # geoscheme = path + "geoscheme.tsv"
    # output = path + 'metadata_geo.tsv'


    focus = ['USA']
    not_focus = ['India', 'South Africa']


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

    # parse subcontinental regions in geoscheme
    scheme_list = open(geoscheme, "r").readlines()[1:]
    geoLevels = {}
    c = 0
    for line in scheme_list:
        if not line.startswith('\n'):
            id = line.split('\t')[2]
            type = line.split('\t')[0]
            if type == 'region_exposure':
                members = line.split('\t')[5].split(',') # elements inside the subarea
                for country in members:
                    iso = get_iso(country.strip())
                    geoLevels[iso] = id

            # parse subnational regions for countries in geoscheme
            if type == 'country_exposure':
                members = line.split('\t')[5].split(',') # elements inside the subarea
                for state in members:
                    if state.strip() not in geoLevels.keys():
                        geoLevels[state.strip()] = id

            # parse subareas for states in geoscheme
            if type == 'division_exposure':
                members = line.split('\t')[5].split(',')  # elements inside the subarea
                for zipcode in members:
                    if zipcode.strip() not in geoLevels.keys():
                        geoLevels[zipcode.strip()] = id
    # print(geoLevels)


    # open metadata file as dataframe
    dfN = pd.read_csv(metadata, encoding='utf-8', sep='\t')
    dfN['region_exposure'] = dfN['iso'].map(geoLevels) # add 'column' region in metadata

    notfound = []
    # convert sets of locations into sub-locations
    print('\nApplying geo-schemes...')
    dfN.fillna('', inplace=True)
    search = SearchEngine(simple_zipcode=True)
    for idx, row in dfN.iterrows():

        # flatten country and division names
        country = dfN.loc[idx, 'country_exposure']
        if country not in focus:
            dfN.loc[idx, 'category'] = 'International'
            dfN.loc[idx, 'division_exposure'] = ''
            dfN.loc[idx, 'location'] = ''


        if country == 'USA':
            dfN.loc[idx, 'category'] = 'USA'

        # convert sets of states into subnational regions
        division = dfN.loc[idx, 'division_exposure']
        category = dfN.loc[idx, 'category']

        # # flatten location names as division names for divisions that are not a focus of study
        # if division not in focus:
        #     dfN.loc[idx, 'location'] = division

        print('Processing metadata for... ' + row['strain'])

    # report errors
    if len(notfound) > 0:
        print('\nSome locations were not assigned to sub-locations, and were not exported. Check for typos.')
        for entry in notfound:
            print('- ' + entry)

    # drop lines if samples not in focal group
    dfN = dfN[~dfN['region_exposure'].isin(not_focus)]
    dfN = dfN[~dfN['country_exposure'].isin(not_focus)]


    # export converted metadata
    dfN.to_csv(output, sep='\t', index=False)
print('\nMetadata file successfully reformatted applying geo-scheme!\n')
