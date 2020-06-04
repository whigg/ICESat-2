import getpass
import os
import socket
import re

###Written by Eduard Heijkoop, University of Colorado###
###Eduard.Heijkoop@colorado.edu###

#This script will generate an NSIDC token, used as authentication to download ICESat-2 data.
#The token will be printed to the screen, stored in "Token.txt" and is valid for 30 days.

user = 'YOUR_NASA_EARTHDATA_USERNAME'
pw = getpass.getpass()
s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
s.connect(("8.8.8.8", 80))
ip_address = s.getsockname()[0]
s.close()


token_command = 'curl -X POST --header \"Content-Type: application/xml\" -d \'<token><username>'+user+'</username><password>'+pw+'</password><client_id>NSIDC_client_id</client_id><user_ip_address>'+ip_address+'</user_ip_address> </token>\' https://cmr.earthdata.nasa.gov/legacy-services/rest/tokens -o tmp_token.txt --silent'

# print(token_command)
os.system(token_command)

for line in open('tmp_token.txt'):
	if "<id>" in line:
		line_req = line

token = re.search('<id>(.*)</id>', line_req)
token = token.group(1)
print(token)


os.system('echo ' + token + ' > Token.txt')
os.system('rm tmp_token.txt')
