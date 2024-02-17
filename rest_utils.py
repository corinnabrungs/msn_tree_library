import requests
import json


def get_json_response(url, post=False, timeout=10):
    """

    :param url:
    :param post:
    :param timeout: default is 10 seconds
    :return:
    """
    try:
        response = (
            requests.post(url, timeout=timeout)
            if post
            else requests.get(url, timeout=timeout)
        )
        response.raise_for_status()
        return json.loads(response.text)
    except requests.exceptions.HTTPError as errh:
        print("Http Error:", errh)
    except requests.exceptions.ConnectionError as errc:
        print("Error Connecting:", errc)
    except requests.exceptions.Timeout as errt:
        print("Timeout Error:", errt)
    except requests.exceptions.RequestException as err:
        print("Other error:", err)
    # on error return None
    return None


def get_json_response_with_headers(url, headers, body, post=False):
    try:
        req_type = "POST" if post else "GET"
        response = requests.request(req_type, url, json=body, headers=headers)
        response.raise_for_status()
        return json.loads(response.text)
    except requests.exceptions.HTTPError as errh:
        print("Http Error:", errh)
    except requests.exceptions.ConnectionError as errc:
        print("Error Connecting:", errc)
    except requests.exceptions.Timeout as errt:
        print("Timeout Error:", errt)
    except requests.exceptions.RequestException as err:
        print("Other error:", err)
    # on error return None
    return None


def json_col(
    result_df, json_column, prefix, field, apply_function=None, new_col_name=None
):
    if new_col_name is None:
        new_col_name = field
    full_column_name = f"{prefix}_{new_col_name}"
    if apply_function is None:
        result_df[full_column_name] = [
            jo[field] if jo and field in jo else None for jo in json_column
        ]
    else:
        result_df[full_column_name] = [
            None if jo is None or field not in jo else apply_function(jo[field])
            for jo in json_column
        ]
    return result_df


def extract_external_descriptors(json_array):
    # [{'source': 'CHEBI', 'source_id': 'CHEBI:48565', 'annotations': ['organic heteropentacyclic compound',
    # 'methyl ester', 'yohimban alkaloid']}]
    return join(
        [
            "{} ({}):{}".format(json["source"], json["source_id"], json["annotations"])
            for json in json_array
        ]
    )


def extract_name(json):
    return json["name"] if json is not None else None


def join(json_array, sep=";"):
    return sep.join(json_array)


def join_by_field(json, field, sep=";"):
    if json is None:
        return None
    return sep.join(json[field])


def join_array_by_field(json_array, field, sep=";"):
    return sep.join([json[field] for json in json_array if json and field in json])


def extract_names_array(json_array):
    return join_array_by_field(json_array, "name")
