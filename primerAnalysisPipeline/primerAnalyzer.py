from argparse import ArgumentParser

from analysisTools import get_primer_combinations, get_access_token, self_dimerization, hetero_dimerization, binding_location_calculator, hairpin_analyzer
import pandas as pd
import time

OFFSET = 15
DELAY = 0.1

#dataframe for dimerization scores
data = {
    "forward_read": [],
    "reverse_read": [],
    "forward_self_dimer": [], # self dimerization scores
    "forward_3_end": [], # analysis for bonds near the 3' end
    "reverse_self_dimer": [],
    "reverse_3_end": [],
    "het_dimer": [],
    "het_3_end": [],
    "combined_score": [], # combined final score to determine pair performance
    "forward_hairpin": [], # score for hairpin sequence analysis
    "reverse_hairpin": []
}

def token_updater(access_token, expire_time):
    if time.time() >= (expire_time - 60):
        print("Access token expired, creating new token")

        new_token, new_token_expiration = get_access_token(idt_identity, idt_secret, username, password)
        print(f"New token created with expiration: {new_token_expiration}")

        return new_token, (time.time() + new_token_expiration)

    return access_token, expire_time

def self_dimer_scoring (primer, combination, token):
    # Find top self-dimerization scores for the chosen combination
    sequence = primer + combination
    self_dimer = self_dimerization(sequence, token)
    top_ends = []
    top_scores = []

    for i in range(top_reads):

        # self-dimerization delta G score for the ith top interaction
        top_self_dimer = self_dimer[i]
        score = top_self_dimer["DeltaG"]

        # for 3' end binding analysis
        binding = top_self_dimer["Bonds"]
        top_padding = top_self_dimer["TopLinePadding"]
        bond_padding = top_self_dimer["BondLinePadding"]
        bottom_padding = top_self_dimer["BottomLinePadding"]
        top_end, bottom_end = binding_location_calculator(
            sequence,
            binding,
            top_padding,
            bond_padding,
            bottom_padding,
            OFFSET
        )
        print(f"self dimerization score for {combination}: {score}")

        # adds ith dimerization 3' end score and dimerization delta G score
        top_ends.append(top_end + bottom_end)
        top_scores.append(score)

    # adds the 3' end and delta-G means of the top n reads to the output
    output_scoring = [sum(top_ends)/top_reads, sum(top_scores)/top_reads, sequence]

    return output_scoring

def hairpin_scoring (primer, combination, token):

    sequence = primer + combination
    hairpin_output = hairpin_analyzer(sequence, token)[0]

    hairpin_score = hairpin_output["deltaG"]

    return hairpin_score

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-f", "--forward", required=True, help="forward binding sequence")
    parser.add_argument("-r", "--reverse", required=True, help="reverse binding sequence")
    parser.add_argument("-a", "--forward_primer", required=True, help="forward primer sequence without binding region" )
    parser.add_argument("-b", "--reverse_primer", required=True, help="reverse primer sequence without binding region" )
    parser.add_argument("-sd", "--self", required=False, default=1, type=int, help="self-dimerization weight")
    parser.add_argument("-hd", "--hetero", required=False, default=1, type=int, help="heterogeneous dimerization weight")
    parser.add_argument("-bd", "--bind", required=False, default=0, type=int, help="binding dimerization location weight")
    parser.add_argument("-hp", "--hairpin", required=False, default=0, type=int, help="hairpin formation weight")
    parser.add_argument("-t", "--topn", type=int, required=False, default=1, help="range of top n binding sequences in dimerization")
    parser.add_argument("-l", "--length", required=True, type=int, help="length of binding sequence")
    parser.add_argument("-i", "--id", required=True, help="id of the IDT account")
    parser.add_argument("-sc", "--secret", required=True, help="secret of the IDT account")
    parser.add_argument("-u", "--username", required=True, help="username of the IDT account")
    parser.add_argument("-p", "--password", required=True, help="password of the IDT account")
    parser.add_argument("-o", "--output", required=True, help="output file")


    args = parser.parse_args()

    # raw reads assignment
    forward_raw = args.forward
    reverse_raw = args.reverse
    forward_primer = args.forward_primer
    reverse_primer = args.reverse_primer
    binding_sequence_length = args.length
    output_file = args.output

    #average top n reads used in the dimerization proces
    top_reads = args.topn

    # accesses IDT API
    idt_identity = args.id
    idt_secret = args.secret
    username = args.username
    password = args.password

    # Creates the access token for your IDT API, with expiry timer
    access_token, token_expiration = get_access_token(idt_identity, idt_secret, username, password)
    expire_time = time.time() + token_expiration
    print(f"access_token created, expires in: {token_expiration}")

    # gets every possible shifted combination of the forward and reverse primers
    forward_combinations = get_primer_combinations(forward_raw, binding_sequence_length)
    reverse_combinations = get_primer_combinations(reverse_raw, binding_sequence_length)

    # weights for the self-dimerization, hetero-dimerization, binding location for dimerization, and hairpin structure formation
    self_weight = args.self
    het_weight = args.hetero
    binding_weight = args.bind
    hairpin_weight = args.hairpin

    dimerization_df = pd.DataFrame(data)

    # dictionary for all forward and reverse dimerization data
    forward_dict = {}
    reverse_dict = {}

    print(f"-----STARTING ON FORWARD READS-----")

    for forward_combination in forward_combinations:
        #checks/updates token if it is about to expire
        access_token, expire_time = token_updater(access_token, expire_time)

        print(f"starting forward: {forward_combination}, time: {time.time()}")

        self_dimer_output = self_dimer_scoring(forward_primer, forward_combination, access_token)
        hairpin_output = hairpin_scoring(forward_primer, forward_combination, access_token)

        forward_output = self_dimer_output + [hairpin_output]

        # updates dictionary with information on new read based on IDT API
        forward_dict[forward_combination] = forward_output

        # Delay for less than 500 API calls per minute
        time.sleep(DELAY)

    print(f"-----STARTING ON REVERSE READS-----")

    for reverse_combination in reverse_combinations:
        # checks/updates access tokens
        access_token, expire_time = token_updater(access_token, expire_time)

        print(f"starting reverse: {reverse_combination}, time: {time.time()}")

        self_dimer_output = self_dimer_scoring(reverse_primer, reverse_combination, access_token)
        hairpin_output = hairpin_scoring(reverse_primer, reverse_combination, access_token)

        reverse_output = self_dimer_output + [hairpin_output]

        # Updates dictionary with new reverse read
        reverse_dict[reverse_combination] = reverse_output

        #Delay for 300 API calls per minute
        time.sleep(DELAY)

    print(f"-----STARTING ON HETERO-DIMERIZATION READS-----")

    # Heterogeneous binding analysis
    for forward_combination in forward_combinations:
        for reverse_combination in reverse_combinations:
            # Checks for access tokens and expiry time
            access_token, expire_time = token_updater(access_token, expire_time)

            print(f"starting hetero: {forward_combination} and {reverse_combination}, time: {time.time()}")

            #Find hetero-dimerization scores for all forward and reverse pairings
            forward_self_end, forward_self_score, forward_self_sequence, forward_hairpin = forward_dict[forward_combination]
            reverse_self_end, reverse_self_score, reverse_self_sequence, reverse_hairpin = reverse_dict[reverse_combination]

            # Uses API for top n heterogeneous binding scores
            het_dimerization = hetero_dimerization(forward_self_sequence, reverse_self_sequence, access_token)
            top_het_ends = []
            top_het_scores = []

            for j in range(top_reads):
                # top ith heterogeneous dimerization interaction
                top_het_dimerization = het_dimerization[j]

                # top ith dimerization delta G scoring
                het_score = top_het_dimerization["DeltaG"]

                # for binding analysis
                het_binding = top_het_dimerization["Bonds"]
                het_top_padding = top_het_dimerization["TopLinePadding"]
                het_bond_padding = top_het_dimerization["BondLinePadding"]
                het_bottom_padding = top_het_dimerization["BottomLinePadding"]
                het_end_forward, het_end_reverse = binding_location_calculator(
                    forward_self_sequence,
                    het_binding,
                    het_top_padding,
                    het_bond_padding,
                    het_bottom_padding,
                    OFFSET
                )

                print(f"hetero dimerization score for {forward_combination} and {reverse_combination}: {het_score}")

                top_het_ends.append(het_end_forward + het_end_reverse)
                top_het_scores.append(het_score)

            #average of the top n scores and ends
            top_het_ends = sum(top_het_ends)/top_reads
            top_het_scores = sum(top_het_scores)/top_reads

            # find final combined scores for each forward and reverse combination.
            combined_score = (self_weight * forward_self_score) + (self_weight * reverse_self_score) + (
                het_weight * top_het_scores)

            # IF bd != 0, factors in whether the dimerization occurs in the last 15 nucleotides or not for both self-dimerization and heterogeneous dimerization for the model
            combined_score -= binding_weight * (self_weight * (forward_self_end + reverse_self_end) +
                                                het_weight * top_het_ends)

            # IF hp != 0, factors in whether hairpin formation in the forward and reverse reads for the model
            combined_score += hairpin_weight * (forward_hairpin + reverse_hairpin)

            # new row in combination
            new_combination = pd.DataFrame([{
                "forward_read": forward_combination,
                "reverse_read": reverse_combination,
                "forward_self_dimer": forward_self_score,
                "forward_3_end": forward_self_end,
                "reverse_self_dimer": reverse_self_score,
                "reverse_3_end": reverse_self_end,
                "het_dimer": top_het_scores,
                "het_3_end": top_het_ends,
                "combined_score": combined_score,
                "forward_hairpin": forward_hairpin,
                "reverse_hairpin": reverse_hairpin,
            }])

            # adds row to the database
            dimerization_df = pd.concat([dimerization_df, new_combination], ignore_index=True)
            #Delay for 300 API calls per minute
            time.sleep(DELAY)

    # Sort data
    print(f"Sorting data")
    sorted_data = dimerization_df.sort_values(by=["combined_score"], ascending=False)

    # convert to .csv format
    print(f"Sorting done. writing to {output_file}")
    sorted_data.to_csv(output_file, index=False)