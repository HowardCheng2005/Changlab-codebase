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

    :param client_id: the client ID for API
    :param client_secret: the client secret for API
    :param idt_username: the username for IDT
    :param idt_password: the password for IDT

    :return: the access token and the amount of time it has before expiry
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
    return body_dict["access_token"], body_dict["expires_in"]

def self_dimerization(sequence, access_token):
    """
    Using the access token and a full primer sequence, returns a list of dictionaries that contain all potential self-dimerizations
    and its scores

    :param sequence: the entire potential primer sequence with illumina adaptor, index, primer binding sequence, heterogeneity spacer, and binding sequence
    :param access_token: the access token for the IDT API
    :return: A list of dictionaries containing each a potential self-dimerization and its scores.
    """
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
    """
    Using an access token and a full primary and secondary sequence, returns all potential heterogeneous dimers and their scores

    :param primary_sequence: top primer sequence containing illumina adaptor, index, primer binding sequence, heterogeneity spacer, and binding sequence
    :param secondary_sequence:  bottom primer sequence containing illumina adaptor, index, primer binding sequence, heterogeneity spacer, and binding sequence
    :param access_token: the access token for the IDT API
    :return: a list of dictionaries containing each a potential heterogeneous dimer and their scores
    """
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

def binding_location_calculator(top_sequence,
                                bond_sequence,
                                top_padding,
                                bond_padding,
                                bottom_padding,
                                offset=15
                                ):
    """
    Calculates the amount of bonds at the 3' end of a self or heterogeneous dimerization based off of the IDT binding format.
    This definition takes in the top and bond sequences from the IDT output, along with their padding, and returns the scores
    for both the top and bottom sequences. Scoring is done by adding full or partial bonds as scores by 2 or 1 respectively.

    :param top_sequence: the top primer sequence within the dimerization reaction
    :param bond_sequence: the bonds that happen between the top and bottom primer sequences within the dimerization
    :param top_padding: the offset/padding of the top sequence to match with the bonds and bottom sequence
    :param bond_padding: the offset/padding of the bonds to match with the top and bottom sequence
    :param bottom_padding: the offset/padding of the bottom sequence to match with the top and bond sequence
    :param offset: the amount of nucleotides from the 3' end to check
    :return: the bond scoring in the 3' end of the top and bottom sequence.
    """

    indexed_bonds = ([0]*bond_padding) + bond_sequence

    top_score = 0
    bottom_score = 0

    # Position of the top 3' end within the bond index
    top_prime_end_pos = top_padding + len(top_sequence) - 1


    # needs to work on indexing
    for i in range(offset):

        # check top 3' end
        top_position = top_prime_end_pos - i
        if 0 <= top_position < len(indexed_bonds):
            top_score += indexed_bonds[top_position]

        # check bottom 3' end
        bottom_position = bottom_padding + i
        if 0 <= bottom_position < len(indexed_bonds):
            bottom_score += indexed_bonds[bottom_position]

    return top_score, bottom_score