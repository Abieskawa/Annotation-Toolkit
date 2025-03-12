import os
import urllib.request  # py3
import urllib.parse
import time
import random
import argparse

homepath = os.getcwd()

def main(resulturl):
    # Parse the input URL to extract "id" and "key"
    parsed_url = urllib.parse.urlparse(resulturl)
    query_params = urllib.parse.parse_qs(parsed_url.query)
    id = query_params.get("id", [None])[0]
    key = query_params.get("key", [None])[0]

    if id and key:
        # Download query.ko file with key parameter appended
        downloadfile = os.path.join(homepath, "query.ko")
        if not os.path.isfile(downloadfile):
            url = f"https://www.genome.jp/tools/kaas/files/dl/{id}/query.ko?key={key}"
            download(downloadfile, url)
        # Download MAP.html file with key parameter appended
        downloadfile = os.path.join(homepath, "MAP.html")
        if not os.path.isfile(downloadfile):
            url = f"https://www.genome.jp/kaas-bin/kaas_main?mode=map&id={id}&key={key}"
            download(downloadfile, url)
        # Download BRITE.html file with key parameter appended
        downloadfile = os.path.join(homepath, "BRITE.html")
        if not os.path.isfile(downloadfile):
            url = f"https://www.genome.jp/kaas-bin/kaas_main?mode=brite&id={id}&key={key}"
            download(downloadfile, url)
    else:
        print("Error: URL must contain both 'id' and 'key' parameters.")
        return

    downloadhierdata()
    downloadmapdata()


def download(downloadfile, url):
    try:
        # Setting a User-Agent header to mimic a browser
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as usock:
            content = usock.read().decode('utf-8')
        with open(downloadfile, 'w') as outputfile:
            outputfile.write(content)
        print(f"Downloaded {url} to {downloadfile}")
    except Exception as e:
        print(f"Error downloading {url}: {e}")


def downloadhierdata():
    keep = 0
    kegpath = os.path.join(homepath, "hier")
    kegghier = os.path.join(homepath, "cleanbrite.html")
    downloadfile = os.path.join(homepath, "BRITE.html")
    if not os.path.isfile(kegghier):
        if os.path.isfile(downloadfile):
            with open(kegghier, 'w') as outputfile, open(downloadfile, 'r') as infile:
                for line in infile:
                    if line.strip():
                        if '/kegg-bin/get_htext?q' in line:
                            keep = 1
                        if '</ul>' in line:
                            keep = 0
                        if keep == 1:
                            outputfile.write(line)
    if not os.path.exists(kegpath):
        os.mkdir(kegpath)

    if os.path.isfile(kegghier):
        counter = 0
        existdata = [x.strip() for x in os.listdir(kegpath) if x.endswith(".keg")]
        with open(kegghier, 'r') as lastfile:
            for line in lastfile:
                if line.strip():
                    parts = line.strip().split('+-p+')
                    if len(parts) > 1:
                        name_part = parts[0].split("/kegg-bin/get_htext?")
                        data_part = parts[1].split('">')
                        if len(name_part) > 1 and (name_part[1] not in existdata):
                            url = f"http://www.genome.jp/kegg-bin/download_htext?htext={name_part[1]}&format=htext&filedir={data_part[0]}"
                            kegpathfile = os.path.join(kegpath, name_part[1])
                            print(f"Downloading to {kegpathfile}")
                            download(kegpathfile, url)
                            counter += 1
                            time.sleep(random.randint(1, 3))
        print(f"updated: {counter}")


def downloadmapdata():
    keep = 0
    mappath = os.path.join(homepath, "map")
    keggmap1 = os.path.join(homepath, "cleanmap.html")
    downloadfile = os.path.join(homepath, "MAP.html")
    if not os.path.isfile(keggmap1):
        if os.path.isfile(downloadfile):
            with open(keggmap1, 'w') as outputfile, open(downloadfile, 'r') as infile:
                for line in infile:
                    if line.strip():
                        if 'show_pathway?ko' in line:
                            keep = 1
                        if 'class="ko_menu"' in line:
                            keep = 0
                        if keep == 1:
                            outputfile.write(line)
    if not os.path.exists(mappath):
        os.mkdir(mappath)

    counter = 0
    mapnum = dict()
    if os.path.isfile(keggmap1):
        with open(keggmap1, 'r') as infile:
            for line in infile:
                if line.strip():
                    parts = line.strip().split('show_pathway?ko')
                    if len(parts) > 1:
                        data = parts[1].split('/')
                        if data:
                            mapnum[data[0].split("@")[0]] = ""
    alldata = sorted(mapnum.keys())
    existdata = [x.strip().replace(".png", "")[3:] for x in os.listdir(mappath) if x.endswith(".png")]
    data = set(alldata) - set(existdata)
    for x in data:
        counter += 1
        print(f"Downloading to {mappath}/map{x}")
        url = f"http://www.genome.jp/kegg-bin/show_pathway?map{x}"
        download(os.path.join(mappath, f"map{x}.html"), url)
        url2 = f"http://www.genome.jp/kegg/pathway/map/map{x}.png"
        try:
            urllib.request.urlretrieve(url2, os.path.join(mappath, f"map{x}.png"))
            print(f"Downloaded image: {url2}")
        except Exception as e:
            print(f"Error downloading image {url2}: {e}")
        time.sleep(random.randint(101,300))
    print(f"updated: {counter}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='downloadKAASresult')
    parser.add_argument('KAASresultURL', type=str, help='KAAS result URL, Usage: python downloadkaasdata.py "URL"')
    args = parser.parse_args()
    resulturl = args.KAASresultURL
    main(resulturl)
    print("Finished.")
