import os
import urllib.request # py3
import time
import random
import argparse

homepath=os.getcwd() 
def main(resulturl):
    data=resulturl.split("id=")
    if len(data)>1:
        id=data[1].split('&')[0].strip() 
        downloadfile=os.path.join(homepath,"query.ko")
        if not os.path.isfile(downloadfile):
            url="https://www.genome.jp/tools/kaas/files/dl/"+id+"/query.ko" 
            download(downloadfile,url)
        downloadfile=os.path.join(homepath,"MAP.html")
        if not os.path.isfile(downloadfile):
            url="https://www.genome.jp/kaas-bin/kaas_main?mode=map&id="+data[1].strip() 
            download(downloadfile,url)
        downloadfile=os.path.join(homepath,"BRITE.html")
        if not os.path.isfile(downloadfile):
            url="https://www.genome.jp/kaas-bin/kaas_main?mode=brite&id="+data[1].strip() 
            download(downloadfile,url)

    downloadhierdata()
    downloadmapdata()
    


def download(downloadfile,url):
        outputfile=open(downloadfile,'w')           
        usock = urllib.request.urlopen(url) # py3
        outputfile.write( usock.read().decode('utf-8')) # py3
        outputfile.flush()
        usock.close()
        outputfile.close()
        #print "Donload:"+ downloadfile

   
def downloadhierdata():
    keep=0
    kegpath=os.path.join(homepath,"hier")
    kegghier=os.path.join(homepath,"cleanbrite.html")
    downloadfile=os.path.join(homepath,"BRITE.html")
    if not os.path.isfile(kegghier):
        if os.path.isfile(downloadfile):
            outputfile=open(kegghier,'w') 
            infile=open(downloadfile,'r') 
            result=[str(x) for x in infile if len(x.strip()) != 0 ]
            for x in result:
                data=x.split('/kegg-bin/get_htext?q')
                if len(data)>1:
                    keep=1
                data=x.split('</ul>')
                if len(data)>1:
                    keep=0
                if keep==1:
                    outputfile.write( x)
            infile.close()
            outputfile.flush()
            outputfile.close()

    if not os.path.exists(kegpath):
        os.mkdir(kegpath)

    if os.path.isfile(kegghier):
        counter=0
        existdata= [str(x).strip() for x in os.listdir(kegpath) if x.endswith(".keg")]
        lastfile=open(kegghier,'r') 
        result=[str(x).strip().split('+-p+') for x in lastfile if len(x.strip()) != 0 ]
        for x in result:
            if len(x)>1:
                name=x[0].split("/kegg-bin/get_htext?")
                data=x[1].split('">') 
                if (str(name[1]) not in existdata):
                    url="http://www.genome.jp/kegg-bin/download_htext?htext="+name[1]+"&format=htext&filedir="+data[0] 
                    kegpathfile=os.path.join(kegpath, str(name[1])) 
                    print(f"Downloading to {kegpathfile}") 
                    download(kegpathfile,url)
                    counter=counter+1
                    time.sleep(random.randint(101,300))
        lastfile.close()
        print(f"updated: {str(counter)}")

def downloadmapdata():
    keep=0
    mappath=os.path.join(homepath,"map")
    keggmap1=os.path.join(homepath,"cleanmap.html")
    downloadfile=os.path.join(homepath,"MAP.html")
    if not os.path.isfile(keggmap1):
        if os.path.isfile(downloadfile):
            outputfile=open(keggmap1,'w') 
            infile=open(downloadfile,'r') 
            result=[str(x) for x in infile if len(x.strip()) != 0 ]
            for x in result:
                data=x.split('show_pathway?ko')
                if len(data)>1:
                    keep=1
                data=x.split('class="ko_menu"')
                if len(data)>1:
                    keep=0
                if keep==1:
                    outputfile.write(x)
            infile.close()
            outputfile.flush()
            outputfile.close()

    if not os.path.exists(mappath):
        os.mkdir(mappath)

    counter=0
    mapnum=dict()
    if os.path.isfile(keggmap1):
        infile=open(keggmap1,'r') 
        result=[str(x).strip().split('show_pathway?ko') for x in infile if len(x.strip()) != 0 ]
        for x in result:
            if len(x)>1:
                data=x[1].split('/')
                mapnum[data[0].split("@")[0]]=""
        infile.close()

    result=sorted(mapnum.keys()) 
    alldata=[x for x in result]
    existdata= [str(x).strip().replace(".png","")[3:] for x in os.listdir(mappath) if x.endswith(".png")]
    data= set(alldata)-set(existdata)

    for x in data:      
        counter=counter+1
        print(f"Downloading to {mappath}/map{x}") # py3
        url = 'http://www.genome.jp/kegg-bin/show_pathway?map'+x
        download(os.path.join(mappath,"map"+str(x)+".html"),url)

        url2 = 'http://www.genome.jp/kegg/pathway/map/map'+x+'.png'
        urllib.request.urlretrieve(url2,os.path.join(mappath,"map"+x+".png")) # py3        
        time.sleep(random.randint(101,300))
        
    print(f"updated: {str(counter)}")

              

if __name__=='__main__':
    parser=argparse.ArgumentParser(description='downloadKAASresult')
    parser.add_argument('KAASresultURL',type=str,help='KAAS result URL ,  Usage:  python downloadkaasdata.py "URL" ')
    args=parser.parse_args()
    resulturl=args.KAASresultURL
    main(resulturl)
    print("Finished.")