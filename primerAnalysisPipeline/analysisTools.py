from base64 import b64encode
import json
from urllib import request, parse

def get_primer_combinations(sequence, length):
    """
    Create a list of possible binding sequence combinations based on a set region and a predetermined length
    :param sequence: the entire potential binding region
    :param length: length of the binding sequence for the primer

    :return: list of binding sequence combinations
    """
    region_length = len(sequence)

    if region_length < length:
        raise ValueError(f"binding region length {region_length} is less than sequence length {length}")

    combination_count = region_length - length + 1
    total_combinations = []

    for i in range(combination_count):

        primer = sequence[i:i+length]
        total_combinations.append(primer)

    print(f"Printed {combination_count} primer combinations")
    return total_combinations


def get_access_token(client_id, client_secret, idt_username, idt_password):
    """
    Create the HTTP request, transmit it, and then parse the response for the
    access token.

    The body_dict will also contain the fields "expires_in" that provides the
    time window the token is valid for (in seconds) and "token_type".
    """

    # Construct the HTTP request
    authorization_string = b64encode(bytes(client_id + ":" + client_secret, "utf-8")).decode()
    request_headers = {"Content-Type": "application/x-www-form-urlencoded",
                       "Authorization": "Basic " + authorization_string}

    data_dict = {"grant_type": "password",
                 "scope": "test",
                 "username": idt_username,
                 "password": idt_password}
    request_data = parse.urlencode(data_dict).encode()

    post_request = request.Request("https://www.idtdna.com/Identityserver/connect/token",
                                   data=request_data,
                                   headers=request_headers,
                                   method="POST")

    # Transmit the HTTP request and get HTTP response
    response = request.urlopen(post_request)

    # Process the HTTP response for the desired data
    body = response.read().decode()

    # Error and return the response from the endpoint if there was a problem
    if response.status != 200:
        raise RuntimeError(f"Request failed with error code: {response.status}. Body: {body}")

    body_dict = json.loads(body)
    return body_dict["access_token"]

def self_dimerization(sequence, access_token):
    url = "https://www.idtdna.com/restapi/v1/OligoAnalyzer/SelfDimer"

    params = parse.urlencode({"primary": sequence})
    full_url = url + "?" + params

    headers = {
        "Authorization": f"Bearer {access_token}",
        "Accept": "application/json"
    }

    req = request.Request(
        url=full_url,
        headers=headers,
        method="POST",
    )

    response = request.urlopen(req)
    body = response.read().decode()

    if response.status != 200:
        raise RuntimeError(f"Request failed with error code: {response.status}. Body: {body}")

    return json.loads(body)

def hetero_dimerization(primary_sequence, secondary_sequence, access_token):
    url = "https://www.idtdna.com/restapi/v1/OligoAnalyzer/HeteroDimer"

    params = parse.urlencode({"primary": primary_sequence, "secondary": secondary_sequence})
    full_url = url + "?" + params

    headers = {
        "Authorization": f"Bearer {access_token}",
        "Accept": "application/json"
    }

    req = request.Request(
        url=full_url,
        headers=headers,
        method="POST",
    )

    response = request.urlopen(req)
    body = response.read().decode()

    if response.status != 200:
        raise RuntimeError(f"Request failed with error code: {response.status}. Body: {body}")

    return json.loads(body)

