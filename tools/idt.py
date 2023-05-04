#!/home/rdkibler/.conda/envs/domesticator_py36/bin/python
import os
from pathlib import Path
import time
import json
import requests
import base64
import argparse
import datetime

from base64 import b64encode
import json
from urllib import request, parse

data_collection_basedir = "/net/shared/idt_dna/sequence_complexity_data_raw"

def vprint(str, verbose=False, **kwargs):
    if verbose:
        print(str, **kwargs)


def use_dir(dir):
    user_info_file = os.path.expanduser(os.path.join(dir, "info.json"))
    # token_file = os.path.expanduser(os.path.join(dir, "token.json"))
    return user_info_file #, token_file


def ask_for_user_data(user_info_file):
    user_info = {}
    idt_url = "https://www.idtdna.com/pages/tools/apidoc"
    msg = (
        f"The first time you use the IDT "
        f"complexity function requires you to "
        f"input your IDT username, password, "
        f"API client ID and API client secret. "
        f"These will be stored securely in "
        f"'{os.path.abspath(user_info_file)}', "
        f"which only you have access to. Before "
        f"you begin, go to {idt_url} and follow "
        f"all the directions under the 'Get "
        f"access to the API' header. The client "
        f"ID and client description can be "
        f"anything. The client secret will be "
        f"generated for you."
    )

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
        with open(user_info_file, "r") as f:
            user_info = json.load(f)

    else:
        # we have to set things up for the first time
        os.makedirs(os.path.dirname(user_info_file), exist_ok=True)

        user_info = ask_for_user_data(user_info_file)

        # do it this weird way just to make sure nobody can see your private stuff
        # simply creates an empty file (like touch)
        Path(user_info_file).touch()
        # make it so only the user can acces this file
        os.chmod(user_info_file, 0o600)
        # now write the secret stuff ;)
        with open(user_info_file, "w+") as f:
            json.dump(user_info, f)

    if "token_file_path" not in user_info:
        user_info["token_file_path"] = os.path.join(
            os.path.dirname(user_info_file), "token.json"
        )

        with open(user_info_file, "w+") as f:
            json.dump(user_info, f)

    #look for the consent_to_data_collection field
    if "consent_to_data_collection" not in user_info:
        print("We want to build a fast IDT score predictor. In order to do")
        print("this we need to collect the DNA sequences you send to IDT ")
        print("and the scores that IDT sends back. This is optional but your")
        print("contribution will help us build a better predictor. We will NOT")
        print("publish or make publically available any of your data. Do you")
        print("consent? (y/n)")

        while True:
            consent = input()
            if consent.lower() == "y":
                user_info["consent_to_data_collection"] = True
                print("Awesome! Thank you for helping us build a better predictor.")
                break
            elif consent.lower() == "n":
                user_info["consent_to_data_collection"] = False
                print("Ok... but if you change your mind you can always change")
                print("this setting in the file:")
                print(os.path.abspath(user_info_file))
                break
            else:
                print("Please enter y or n")

        with open(user_info_file, "w+") as f:
            json.dump(user_info, f)

    return user_info


def delete_stored_token(token_file):
    if os.path.exists(token_file):
        os.remove(token_file)


def store_token(token, token_file):
    delete_stored_token(token_file)
    # do it this weird way just to make sure nobody can see your private stuff
    # creates an empty file
    Path(token_file).touch()
    # make it so only the user can acces this file
    os.chmod(token_file, 0o600)
    # now write the secret stuff ;)
    with open(token_file, "w+") as f:
        json.dump(token, f)


# IDT's code, because my code is bad
def get_new_token(user_info, verbose=False):
    """
    Create the HTTP request, transmit it, and then parse the response for the
    access token.

    The body_dict will also contain the fields "expires_in" that provides the
    time window the token is valid for (in seconds) and "token_type".
    """
    vprint("getting new token", verbose)
    client_id = user_info["ID"]
    client_secret = user_info["secret"]
    idt_username = user_info["username"]
    idt_password = user_info["password"]

    # Construct the HTTP request
    authorization_string = b64encode(
        bytes(client_id + ":" + client_secret, "utf-8")
    ).decode()
    request_headers = {
        "Content-Type": "application/x-www-form-urlencoded",
        "Authorization": "Basic " + authorization_string,
    }

    data_dict = {
        "grant_type": "password",
        "scope": "test",
        "username": idt_username,
        "password": idt_password,
    }
    request_data = parse.urlencode(data_dict).encode()

    post_request = request.Request(
        "https://www.idtdna.com/Identityserver/connect/token",
        data=request_data,
        headers=request_headers,
        method="POST",
    )

    # Transmit the HTTP request and get HTTP response
    response = request.urlopen(post_request)

    # Process the HTTP response for the desired data
    body = response.read().decode()

    # Error and return the response from the endpoint if there was a problem
    if response.status != 200:
        raise RuntimeError(
            "Request failed with error code:" + response.status + "\nBody:\n" + body
        )

    body_dict = json.loads(body)
    return body_dict


def get_stored_token(token_file):
    with open(token_file, "r") as f:
        return json.load(f)


def get_token(user_info, verbose=False):
    token_file = user_info['token_file_path']
    get_token_flag = False
    if os.path.exists(token_file):
        vprint(f"using token stored at {token_file}", verbose)
        modified_time = os.path.getmtime(token_file)
        current_time = time.time()
        token = get_stored_token(token_file)
        if current_time - modified_time > token["expires_in"]:
            vprint("token expired", verbose)
            get_token_flag = True
    else:
        vprint(f"no file found at {token_file}", verbose)
        get_token_flag = True

    if get_token_flag:
        token = get_new_token(user_info, verbose)
        vprint(f"storing token at {token_file}", verbose)
        store_token(token, token_file)
    return token

def store_response(response_dict, seq, kind):
    date = datetime.datetime.now().strftime("%Y-%m-%d")
    timestamp = datetime.datetime.now().strftime("%H-%M-%S-%f")[:-3]
    username = os.getlogin()

    data_collection_dir = os.path.join(data_collection_basedir,username, date)
    os.makedirs(data_collection_dir)

    data = {}
    data['response'] = response_dict
    data['timestamp'] = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    data['username'] = username
    data['kind'] = str(kind)
    data['sequence'] = str(seq)

    filename = f"{timestamp}.json"
    filepath = os.path.join(data_collection_dir, filename)

    with open(filepath, "w+") as f:
        json.dump(data, f)



def query_complexity(seq, user_info, verbose=False, kind='gene'):
    token = get_token(user_info)['access_token']

    url_dict = {}
    url_dict['gene'] = "https://www.idtdna.com/Restapi/v1/Complexities/ScreenGeneSequences"
    url_dict['gblock'] = "https://www.idtdna.com/Restapi/v1/Complexities/ScreenGblockSequences"
    url_dict['gblock_hifi'] = "https://www.idtdna.com/Restapi/v1/Complexities/ScreenGblockHifiSequences"
    url_dict['eblock'] = "https://www.idtdna.com/Restapi/v1/Complexities/ScreenEblockSequences"
    url_dict['old'] = "https://www.idtdna.com/api/complexities/screengBlockSequences"

    url = url_dict[kind]

    vprint(f"querying complexity", verbose)

    payload = f'[{{"Name":"My gBlock","Sequence":"{seq}"}}]'
    headers = {"Content-Type": "application/json", "Authorization": f"Bearer {token}"}

    while True:
        response = requests.request("POST", url, headers=headers, data=payload)
        if response.status_code == 429:
            print(response.text)
            vprint("Rate limited, waiting 10 seconds", verbose)
            time.sleep(10)
        else:
            break

    if int(response.status_code) != 200:

        print(response)
        print(dir(response))
        print("JSON")
        print(response.json())
        print("TEXT")
        print(response.text)
        print("CONTENT")
        print(response.content)
        print("STATUS_CODE", response.status_code)
        print("HEADERS")
        print(response.headers)
        print("COOKIES")
        print(response.cookies)
        print("URL")
        print(response.url)
        print("HISTORY")
        print(response.history)
        print("ENCODING")
        print(response.encoding)
        print("REASON")
        print(response.reason)
        print("OK")

        raise RuntimeError(
            f"Request failed with error code: {response.status_code} \nBody:\n{response.text}"
        )
    
    response_dict = json.loads(response.text)
    print(response_dict)

    if user_info['consent_to_data_collection']:
        vprint(f"storing response at {data_collection_basedir}", verbose)
        store_response(response_dict, seq, kind)

    return response_dict


if __name__ == "__main__":
    import argparse
    from Bio import SeqIO

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "fasta",
        type=str,
        help="A fasta file containing the sequences you want to check for complexity",
    )
    parser.add_argument("--credential_dir", default="~/idt_credentials")
    parser.add_argument('-k', '--kind', type=str, help='kind of sequence to query', default='gene', choices = ['gene','gblock','gblock_hifi','eblock','old'])
    args = parser.parse_args()

    user_info_file = use_dir(args.credential_dir)
    idt_user_info = get_user_info(user_info_file)

    for record in SeqIO.parse(args.fasta, "fasta"):
        response = query_complexity(record.seq, idt_user_info, kind=args.kind)[0]
        score = 0
        if len(response) == 0:
            score = 0
        else:
            for issue in response:
                score += issue["Score"]
        
        print(record.id, record.seq, score)

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
