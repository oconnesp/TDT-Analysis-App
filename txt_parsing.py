import re
import numpy as np

def extract_from_txt (filename, no_trials):
    responses = []
    ISIs = []

    response_pattern = r"Response:\s*(\d(?:,\s*\d)*)"
    ISI_pattern = r"TDTArray:\s*([0-9eE.+-]+(?:,\s*[0-9eE.+-]+)*)"

    with open (filename, "r") as file:
        for line in file:
            resp_match = re.search (response_pattern, line)
            ISI_match = re.search (ISI_pattern, line)
            if resp_match:
                resp_nums = np.array(list(map(int, resp_match.group(1).split(','))))
                responses.append(resp_nums)
            if ISI_match:
                ISI_nums = np.array([float(n) for n in ISI_match.group(1).split(',')]) * 1000
                ISIs.append((ISI_nums))
                if len(ISIs) == no_trials:
                    break
    all_isis = np.unique(np.concatenate(ISIs))
    return ISIs, responses, all_isis
