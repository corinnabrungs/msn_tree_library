def map_clinical_phase_to_number(phase):
    match (str(phase)):
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
        case "4" | "4.0" | "Phase 4" | "Launched":
            return 4
        case _:
            return phase


def get_clinical_phase_description(number):
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
