def combine_polarity(old, new):
    if old == "both":
        return "both"
    match new:
        case "both":
            return new
        case "positive":
            return "both" if old == "negative" else "positive"
        case "negative":
            return "both" if old == "positive" else "negative"
        case _:
            return old
