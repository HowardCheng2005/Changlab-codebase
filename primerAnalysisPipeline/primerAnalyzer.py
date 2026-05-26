from argparse import ArgumentParser
from analysisTools import get_primer_combinations, get_access_token, self_dimerization, hetero_dimerization
import pandas as pd
import time

OFFSET = 15

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-f", "--forward", required=True, help="forward binding sequence")
    parser.add_argument("-r", "--reverse", required=True, help="reverse binding sequence")
    parser.add_argument("-a", "--forward_primer", required=True, help="forward primer sequence without binding region" )
    parser.add_argument("-b", "--reverse_primer", required=True, help="reverse primer sequence without binding region" )
    parser.add_argument("-sd", "--self", required=False, default=1, type=int, help="self-dimerization weight")
    parser.add_argument("-hd", "--hetero", required=False, default=1, type=int, help="heterogeneous dimerization weight")
    parser.add_argument("-bd", "--bind", required=False, default=0, type=int, help="binding dimerization location weight")
    parser.add_argument("-l", "--length", required=True, type=int, help="length of binding sequence")
    parser.add_argument("-i", "--id", required=True, help="id of the IDT account")
    parser.add_argument("-sc", "--secret", required=True, help="secret of the IDT account")
    parser.add_argument("-u", "--username", required=True, help="username of the IDT account")
    parser.add_argument("-p", "--password", required=True, help="password of the IDT account")
    parser.add_argument("-o", "--output", required=True, help="output file")


    args = parser.parse_args()

    forward_sequence = args.forward
    reverse_sequence = args.reverse
    forward_primer = args.forward_primer
    reverse_primer = args.reverse_primer
    binding_sequence_length = args.length
    output_file = args.output

    # accesses IDT API
    idt_identity = args.id
    idt_secret = args.secret
    username = args.username
    password = args.password

    # access token for IDT
    access_token = get_access_token(idt_identity, idt_secret, username, password)

    # gets every possible shifted combination of the forward and reverse primers
    forward_combinations = get_primer_combinations(forward_sequence, binding_sequence_length)
    reverse_combinations = get_primer_combinations(reverse_sequence, binding_sequence_length)

    # weights for the self-dimerization, hetero-dimerization, and binding location for dimerization
    self_weight = args.self
    het_weight = args.hetero
    binding_weight = args.bind

    #dataframe for dimerization scores
    data = {
        "forward_read": [],
        "reverse_read": [],
        "forward_self_dimer": [],
        "forward_binding_sequence": [],
        "reverse_self_dimer": [],
        "reverse_binding_sequence": [],
        "het_dimer": [],
        "het_binding_sequence": [],
        "combined_score": []
    }
    dimerization_df = pd.DataFrame(data)

    # dictionary for all forward and reverse dimerization data
    forward_dict = {}
    reverse_dict = {}

    for forward_combination in forward_combinations:
        # Find self-dimerization scores for all forward combinations
        forward_sequence = forward_primer + forward_combination
        forward_self_dimer = self_dimerization(forward_sequence, access_token)[0]
        forward_score = forward_self_dimer["DeltaG"]
        forward_binding_sequence = forward_self_dimer["Bonds"]
        print(f"Forward self dimerization score for {forward_combination}: {forward_score}")

        forward_dict[forward_combination] = [forward_binding_sequence, forward_score, forward_sequence]
        # Delay for 300 API calls per minute
        time.sleep(0.2)


    for reverse_combination in reverse_combinations:
        # Find self-dimerization scores for all reverse combinations
        reverse_sequence = reverse_primer + reverse_combination
        reverse_self_dimerization = self_dimerization(reverse_sequence, access_token)[0]
        reverse_score = reverse_self_dimerization["DeltaG"]
        reverse_binding_sequence = reverse_self_dimerization["Bonds"]
        print(f"Reverse self dimerization score for {reverse_combination}: {reverse_score}")

        reverse_dict[reverse_combination] = [reverse_binding_sequence, reverse_score, reverse_sequence]
        #Delay for 300 API calls per minute
        time.sleep(0.2)


    for forward_combination in forward_combinations:
        for reverse_combination in reverse_combinations:
            #Find hetero-dimerization scores for all forward and reverse pairings
            forward_self_binding, forward_self_score, forward_self_sequence = forward_dict[forward_combination]
            reverse_self_binding, reverse_self_score, reverse_self_sequence = reverse_dict[reverse_combination]

            het_dimerization = hetero_dimerization(forward_self_sequence, reverse_self_sequence, access_token)[0]
            het_score = hetero_dimerization["DeltaG"]
            het_binding_sequence = hetero_dimerization["Bonds"]

            print(f"hetero dimerization score for {forward_combination} and {reverse_combination}: {het_score}")

            # find final combined scores for each forward and reverse combination.
            combined_score = (self_weight * forward_self_score) + (self_weight * reverse_self_score) + (
                het_weight * het_score)

            # IF bd != 0, factors in whether the dimerization occurs in the last 15 nucleotides or not for both self-dimerization and heterogeneous dimerization


            new_combination = pd.DataFrame([{
                "forward_read": forward_combination,
                "reverse_read": reverse_combination,
                "forward_self_dimer": forward_self_score,
                "forward_binding_sequence": forward_self_binding,
                "reverse_self_dimer": reverse_self_score,
                "reverse_binding_sequence": reverse_self_binding,
                "het_dimer": het_score,
                "het_binding_sequence": het_binding_sequence,
                "combined_score": combined_score
            }])

            dimerization_df = pd.concat([data, new_combination], ignore_index=True)
            #Delay for 300 API calls per minute
            time.sleep(0.2)

    # Sort data
    print(f"Sorting data")
    sorted_data = dimerization_df.sort_values(by=["combined_score"], ascending=False)

    # convert to .csv format
    print(f"Sorting done. writing to {output_file}")
    sorted_data.to_csv(output_file, index=False)