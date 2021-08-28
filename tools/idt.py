#!/home/rdkibler/.conda/envs/domesticator_py36/bin/python
import os
from pathlib import Path
import time
import json
import requests
import base64
import argparse

idt_dir="~/.domesticator"

user_info_file = os.path.expanduser(os.path.join(idt_dir, "info.json"))
token_file = os.path.expanduser(os.path.join(idt_dir, "token.json"))


def vprint(str,verbose=False,**kwargs):
	if verbose: print(str,**kwargs)

def use_dir(dir):
	user_info_file = os.path.expanduser(os.path.join(dir, "info.json"))
	token_file = os.path.expanduser(os.path.join(dir, "token.json"))
	return user_info_file, token_file

def ask_for_user_data(user_info_file):
	user_info = {}
	idt_url = "https://www.idtdna.com/pages/tools/apidoc"
	msg = f"The first time you use the IDT " \
		f"complexity function requires you to " \
		f"input your IDT username, password, " \
		f"API client ID and API client secret. "\
		f"These will be stored securely in " \
		f"'{os.path.abspath(user_info_file)}', " \
		f"which only you have access to. Before " \
		f"you begin, go to {idt_url} and follow " \
		f"all the directions under the 'Get " \
		f"access to the API' header. The client " \
		f"ID and client description can be " \
		f"anything. The client secret will be " \
		f"generated for you."

	print(msg)
	print("1) Please enter your IDT account username: ")
	user_info["username"] = input()
	print("2) Please enter your IDT account password: ")
	user_info["password"] = input()
	print("3) Please enter you API client ID: ")
	user_info["ID"] = input()
	print("4) Please enter your API secret: ")
	user_info["secret"] = input()

	return user_info

def get_user_info(user_info_file):
	if os.path.exists(user_info_file):
		with open(user_info_file,'r') as f:
			user_info = json.load(f)

	else:
		#we have to set things up for the first time
		os.makedirs(os.path.dirname(user_info_file),exist_ok=True)

		user_info = ask_for_user_data(user_info_file)

		#do it this weird way just to make sure nobody can see your private stuff
		#simply creates an empty file (like touch)
		Path(user_info_file).touch()
		#make it so only the user can acces this file
		os.chmod(user_info_file,0o600)
		#now write the secret stuff ;)
		with open(user_info_file,'w+') as f:
			json.dump(user_info,f)

	return user_info

def delete_stored_token(token_file):
	if os.path.exists(token_file):
		os.remove(token_file)

def store_token(token,token_file):
	delete_stored_token(token_file)
	#do it this weird way just to make sure nobody can see your private stuff
	#creates an empty file
	Path(user_info_file).touch()
	#make it so only the user can acces this file
	os.chmod(user_info_file,0o600)
	#now write the secret stuff ;)
	with open(token_file,'w+') as f:
		json.dump(token,f)


def get_new_token(user_info, verbose=False):
	vprint("getting new token", verbose)
	client_info = user_info["ID"]+":"+user_info["secret"]
	byte_encoded_client_info = client_info.encode("utf-8")
	base64_client_info = base64.urlsafe_b64encode(byte_encoded_client_info).decode('utf8')

	url = "https://www.idtdna.com/Identityserver/connect/token"

	payload = f'grant_type=password&username={user_info["username"]}&password={user_info["password"]}&scope=test'
	headers = {
		'Content-Type': 'application/x-www-form-urlencoded',
		'Authorization': f'Basic {base64_client_info}'
	}
	tries = 5
	while(tries > 0):
		response = requests.request("POST", url, headers=headers, data = payload)

		response_dict = json.loads(response.text)
		tries =- 1
		if "access_token" in response_dict:
			vprint("token acquired",verbose)
			break
		else:
			for key in response_dict.keys():
				print(key, ":", response_dict[key])
			vprint(f"ERROR: {response_dict['Message']}\nTrying again in 5 seconds...",verbose)
			time.sleep(5)
	if "access_token" not in response_dict:
		raise RuntimeError("Could not get access token")
	return response_dict

def get_stored_token(token_file):
	with open(token_file,'r') as f:
		return json.load(f)

def get_token(token_file, user_info,verbose=False):
	get_token_flag = False
	if os.path.exists(token_file):
                vprint(f"using token stored at {token_file}",verbose)
                modified_time = os.path.getmtime(token_file)
                current_time = time.time()
                token = get_stored_token(token_file)
                if current_time - modified_time > token["expires_in"]:
                        vprint("token expired",verbose)
                        get_token_flag = True
	else:
                vprint(f"no file found at {token_file}",verbose)
                get_token_flag = True

	if get_token_flag:
                token = get_new_token(user_info,verbose)
                vprint(f"storing token at {token_file}",verbose)
                store_token(token,token_file)
	return token

def query_complexity(seq, token,verbose=False):
	vprint("querying complexity",verbose)
	url = "https://www.idtdna.com/api/complexities/screengBlockSequences"

	payload = f'[{{"Name":"My gBlock","Sequence":"{seq}"}}]'
	headers = {
		'Content-Type': 'application/json',
		'Authorization': f'Bearer {token}'
	}

	response = requests.request("POST", url, headers=headers, data = payload)

	return json.loads(response.text)

if __name__ == "__main__":
	import argparse
	from Bio import SeqIO

	parser = argparse.ArgumentParser()
	parser.add_argument("fasta", type=str, help="A fasta file containing the sequences you want to check for complexity")
	parser.add_argument("--credential_dir", default="~/.domesticator")
	args = parser.parse_args()

	user_info_file, token_file = use_dir(args.credential_dir)
	idt_user_info = get_user_info(user_info_file)


	for record in SeqIO.parse(args.fasta,"fasta"):
		token = get_token(token_file, idt_user_info)
		response = query_complexity(record.seq, token["access_token"])[0]
		score = 0
		if len(response) == 0:
			score = 0
		else:
			for issue in response:
				score += issue["Score"]
		print(record.id,record.seq,score)

"""
Accepted - Moderate Complexity (Scores between 7 and 19)

Some complexities exist that may interfere with or delay manufacturing. If it is possible to reduce these complexities please do so, otherwise we will attempt this order.
"""




"""
ForwardLocations is a list of problem locations fwd. I believe these are zero index
ReverseLocations is a list of problem locations rev. I believe these are zero index
RepeatedSegment has the problem sequence (just need the length from this)

things with lists:
SSA Repeat 3
SSA Repeat 5
Pseudo Terminal Repeat
Single Repeat Percentage
Single Repeat Overall Bases
Repeat Length



things with other stuff:
Overall Repeat (might have to parse msg?)


StartIndex has the start (duH)
Not sure how to determine stopping point unless parse msg

things with index
Windowed Repeat Percentage

"""
