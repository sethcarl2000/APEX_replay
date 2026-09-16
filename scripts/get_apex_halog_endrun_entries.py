
import json
import requests

BASE_URL = "https://logbooks.jlab.org/api/elog/entries"

def query_entries(title_pattern, params: dict, output_fields):

    params["field"] = output_fields
        
    resp = requests.get(BASE_URL, params=params, timeout=30)
    print("Requested:", resp.url)  # handy for checking the encoded URL
    resp.raise_for_status()        # raises an exception on 4xx/5xx errors
    return resp.json()             # parses the JSON into dicts/lists

if __name__ == "__main__":

    fields = ["lognumber", "title", "author", "tags", "created", "body"]
    title_pattern = "End_of_Run_"

    tags = ["EndOfRun", "Autolog"]

    start_date = "2019-02-13T16:19"
    end_date = "2019-03-19T05:16" 
    
    params = {
        "title": "End_of_Run_",
        "tags": tags,
        "startdate": start_date,
        "enddate": end_date,
        "limit": 2000,
        "author": "adaq" 
    }

    data = query_entries(title_pattern, params, fields)
    print(json.dumps(data, indent=2))
