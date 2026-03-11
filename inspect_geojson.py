import urllib.request, json

url = 'https://raw.githubusercontent.com/datasets/geo-boundaries-world-110m/master/countries.geojson'
with urllib.request.urlopen(url, timeout=10) as f:
    gj = json.load(f)

for feat in gj['features']:
    props = feat.get('properties', {})
    if 'philippines' in str(props).lower():
        print('FOUND', props)
        print('type', feat['geometry']['type'])
        break
else:
    print('not found')
    for i, feat in enumerate(gj['features'][:10]):
        print(i, feat['properties'].keys())
