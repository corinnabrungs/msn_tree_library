def combine_smarts(smarts):
    combined = ",".join([f"$({sm})" for sm in smarts])
    return f"[{combined}]"
