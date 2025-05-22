#!/usr/bin/env python3
import os
import urllib.request  # py3
import urllib.parse
import time
import random
import argparse
import re  # for parsing the cleanbrite links
import socket
from urllib.error import URLError

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

    downloadhierdata(id, key)  # Pass id and key parameters
    downloadmapdata()


def download(downloadfile, url, max_retries=5, initial_delay=2):
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
        'Accept-Language': 'en-US,en;q=0.5',
        'Connection': 'keep-alive',
    }
    
    delay = initial_delay
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=30) as usock:
                content = usock.read().decode('utf-8')
            with open(downloadfile, 'w') as outputfile:
                outputfile.write(content)
            print(f"Downloaded {url} to {downloadfile}")
            return True
        except (URLError, socket.error) as e:
            if "Connection reset by peer" in str(e):
                print(f"Connection reset on attempt {attempt+1}/{max_retries} for {url}: {e}")
            else:
                print(f"Error on attempt {attempt+1}/{max_retries} for {url}: {e}")
            
            if attempt < max_retries - 1:
                # Add jitter to the delay
                actual_delay = delay * (0.5 + random.random())
                print(f"Retrying in {actual_delay:.2f} seconds...")
                time.sleep(actual_delay)
                # Exponential backoff
                delay *= 2
            else:
                print(f"Failed to download {url} after {max_retries} attempts")
                return False
        except Exception as e:
            print(f"Unexpected error downloading {url}: {e}")
            if attempt < max_retries - 1:
                time.sleep(delay)
                delay *= 2
            else:
                return False


def download_image(url, filepath, max_retries=5, initial_delay=2):
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
    }
    
    delay = initial_delay
    for attempt in range(max_retries):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=30) as response:
                with open(filepath, 'wb') as out_file:
                    out_file.write(response.read())
            print(f"Downloaded image: {url}")
            return True
        except (URLError, socket.error) as e:
            if attempt < max_retries - 1:
                actual_delay = delay * (0.5 + random.random())
                print(f"Image download error on attempt {attempt+1}/{max_retries}: {e}")
                print(f"Retrying in {actual_delay:.2f} seconds...")
                time.sleep(actual_delay)
                delay *= 2
            else:
                print(f"Failed to download image {url} after {max_retries} attempts: {e}")
                return False
        except Exception as e:
            print(f"Unexpected error downloading image {url}: {e}")
            if attempt < max_retries - 1:
                time.sleep(delay)
                delay *= 2
            else:
                return False


def downloadhierdata(id, key):  # Add id and key parameters
    kegpath     = os.path.join(homepath, "hier")
    raw_brite   = os.path.join(homepath, "BRITE.html")
    clean_brite = os.path.join(homepath, "cleanbrite.html")

    # 1) Extract only the get_htext lines
    if not os.path.isfile(clean_brite) and os.path.isfile(raw_brite):
        with open(raw_brite, 'r') as inp, open(clean_brite, 'w') as outp:
            keep = False
            for line in inp:
                if '/kegg-bin/get_htext?' in line:
                    keep = True
                if keep:
                    outp.write(line)
                if keep and '</ul>' in line:
                    keep = False

    # 2) Ensure the output directory exists
    os.makedirs(kegpath, exist_ok=True)

    # 3) Regex to capture the keg file name only
    pattern = re.compile(r'get_htext\?([^+"]+)')  # Extract just the file name
    existing = {f for f in os.listdir(kegpath) if f.endswith(".keg")}

    counter = 0
    with open(clean_brite, 'r') as brite:
        for line in brite:
            m = pattern.search(line)
            if not m:
                continue

            htext = m.group(1)  # e.g. "q00001.keg"
            if htext in existing:
                continue

            # Updated URL format with id and key parameters
            url = f"https://www.genome.jp/tools/kaas/files/log/result/{id}/{htext}?key={key}"
            outp = os.path.join(kegpath, htext)
            print(f"Downloading to {outp}")
            download(outp, url)
            counter += 1
            time.sleep(random.randint(1, 3))

    print(f"updated: {counter}")


def downloadmapdata():
    keep = 0
    mappath     = os.path.join(homepath, "map")
    keggmap1    = os.path.join(homepath, "cleanmap.html")
    downloadfile = os.path.join(homepath, "MAP.html")

    if not os.path.isfile(keggmap1) and os.path.isfile(downloadfile):
        with open(keggmap1, 'w') as outputfile, open(downloadfile, 'r') as infile:
            for line in infile:
                if line.strip():
                    if 'show_pathway?ko' in line:
                        keep = 1
                    if 'class="ko_menu"' in line:
                        keep = 0
                    if keep == 1:
                        outputfile.write(line)

    os.makedirs(mappath, exist_ok=True)

    counter = 0
    mapnum  = {}
    if os.path.isfile(keggmap1):
        with open(keggmap1, 'r') as infile:
            for line in infile:
                if line.strip():
                    parts = line.strip().split('show_pathway?ko')
                    if len(parts) > 1:
                        data = parts[1].split('/')
                        if data:
                            mapnum[data[0].split("@")[0]] = ""

    alldata   = sorted(mapnum.keys())
    existdata = [x.replace(".png","")[3:] for x in os.listdir(mappath) if x.endswith(".png")]
    data      = set(alldata) - set(existdata)

    for x in data:
        counter += 1
        print(f"Downloading to {mappath}/map{x}")
        url = f"http://www.genome.jp/kegg-bin/show_pathway?map{x}"
        download(os.path.join(mappath, f"map{x}.html"), url)
        url2 = f"http://www.genome.jp/kegg/pathway/map/map{x}.png"
        download_image(url2, os.path.join(mappath, f"map{x}.png"))
        time.sleep(random.randint(1,3))

    print(f"updated: {counter}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='downloadKAASresult')
    parser.add_argument('KAASresultURL', type=str,
                        help='KAAS result URL, Usage: python downloadkaasdata.py "URL"')
    args = parser.parse_args()
    resulturl = args.KAASresultURL
    main(resulturl)
    print("Finished.")