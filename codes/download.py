from bs4 import BeautifulSoup
import urllib.request

data_url = 'https://filesender.surf.nl/?s=download&token=c0d5a284-0754-45bc-9031-ea65d60e96ed'
html = urllib.request.urlopen(data_url)
soup = BeautifulSoup(html.read())
for fq_file in soup.findAll("div", {"class": "file"}):
    fn = fq_file.attrs['data-name']
    link = fq_file.find("a", {'class':'download'}).attrs['href']
    print('curl -o %s %s' %(fn, link))
