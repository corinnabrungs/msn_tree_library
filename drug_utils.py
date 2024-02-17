from pandas_utils import isnull, isnull_or_empty


def map_clinical_phase_to_number(phase):
    if isnull_or_empty(phase):
        return 0

    cleanphase = str(phase).split("/")[-1].strip()
    match (cleanphase):
        case "" | "None" | "NaN" | "nan" | "np.nan" | "0" | "No Development Reported":
            return 0
        case "Preclinical":
            return 0.5
        case "1" | "1.0" | "Phase 1":
            return 1
        case "2" | "2.0" | "Phase 2":
            return 2
        case "3" | "3.0" | "Phase 3":
            return 3
        case "4" | "4.0" | "Phase 4" | "Launched" | "Withdrawn":
            return 4
        case _:
            return 0


def get_clinical_phase_description(number):
    if isnull_or_empty(number):
        return number

    match (str(number)):
        case "" | "None" | "NaN" | "nan" | "np.nan" | "0" | "0.0":
            return ""
        case "0.5":
            return "Preclinical"
        case "1.0" | "1":
            return "Phase 1"
        case "2.0" | "2":
            return "Phase 2"
        case "3.0" | "3":
            return "Phase 3"
        case "4.0" | "4":
            return "Launched"
        case _:
            return number
