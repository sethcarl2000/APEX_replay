import csv
import json
import argparse
import re

BODY_PATTENS: dict[str, re.Pattern] = {
    "run_number" : re.compile(r"Run\s*Number:\s*(\d{4})"),
    "timestamp" : re.compile(r"<h4>Log\s*entry\s*time\s*(\d{2}:\d{2}:\d{2})\s*on\s*(February|March)\s*(\d{1,2}),\s*(\d{4})"),
    "run_type" : re.compile(r"Run_type=(\w+)"),
    "target_type" : re.compile(r"target_type=([a-zA-Z0-9 .%]+),"),
    "comment" : re.compile(r"comment_text=([a-zA-Z0-9 .%;=,]+)"),
    "deadtime_fraction" : re.compile(r"DEAD\s*TIME:\s*([0-9.]+)%"),
    "total_events" : re.compile(r"EVENTS\s*:\s*(\d+)"),
    "run_time_seconds" : re.compile(r"TIME\s*:\s*([0-9.]+)\s*mins"),
    "beam_energy_MeV" : re.compile(r"Tiefenbach\s*6GeV\s*Beam\s*energy\s*\(MeV\)\s*:\s*([0-9.]+)"),
    "beam_current_uA" : re.compile(r"Beam\s*Current\s*:\s*([0-9.]+)"),
    "momentum_RHRS_MeV" : re.compile(r"Right\s*arm\s*momentum\s*:\s*([0-9.])"),
    "momentum_LHRS_MeV" : re.compile(r"Left\s*arm\s*momentum\s*:\s*([0-9.])"),
}



def parse_special_rule(
        name: str, match: re.Match 
) -> str | None:
    
    if name == "timestamp": 
        """Takes a successful 'timestamp' re match and returns a YYYY-MM-DDTHH:MM:SS format"""
            
        # thhe HH:MM:SS should already be formatted correctly
        time = match.group(1).strip()

        month = match.group(2).strip()
        mm = 0


        if month == "February":
            mm = 2
        elif month == "March": 
            mm = 3 
        else: 
            raise ValueError(f"Illegal month name: {month}")

        dd   = int(match.group(3).strip())
        yyyy = int(match.group(4).strip())

        # return the formatted time, making sure to add leading zeros where necessary 
        return f"{yyyy:4d}-{mm:02d}-{dd:02d}T{time}"

    elif name == "run_time_seconds": 
        """Takes elapsed run time, in minutes, and converts it to seconds"""
        seconds = float(match.group(1).strip()) * 60
        return f"{seconds:.2f}"

    elif name == "deadtime_fraction": 
        """Convert this from percent to a fraction"""
        val = float(match.group(1).strip()) / 100. 
        return f"{val:.5f}"

    elif name == "momentum_RHRS_MeV": 
        """Convert this from GeV to MeV"""
        val = float(match.group(1).strip()) * 1000. 
        return f"{val:.2f}"

    elif name == "momentum_LHRS_MeV": 
            """Convert this from GeV to MeV"""
            val = float(match.group(1).strip()) * 1000. 
            return f"{val:.2f}"

    else: 
        """this format is not a special format"""
        return None

    

def extract_fields(
        text: str, 
        patterns: dict[str, re.Pattern]
) -> dict[str, str | None]: 
    """Search for matches of 'patterns' in 'text'"""

    results = {}

    for name, regex in patterns.items(): 
        match = regex.search(text)

        if match: 
            if (special_format := parse_special_rule(name,match)) is not None: 
                # handle the rules which need special formatting
                results[name] = special_format
            else:
                # no special formatting rule here. 
                results[name] = match.group(1).strip()
        else:
            results[name] = None 

    # add some extra results

    # compute how many uC (micro-coulombs)
    if results["beam_current_uA"] is not None and results["run_time_seconds"] is not None: 
        current_uA = float(results["beam_current_uA"])
        time_s = float(results["run_time_seconds"])

        # results["accumulated_charge_uC"] = f"{current_uA * time_s:.4f}" 

    return results


def json_to_csv(
    input_path: str, 
    output_path: str, 
    patterns: dict[str, re.Pattern] = BODY_PATTENS, 
    require_all: bool = True
) -> tuple[int, int]: 
    """Take HALOG input json, and return the desired data fields into a csv
    
    Returns(rows_written, rows_skipped)
    """

    written =0 
    skipped =0 
    
    print(f"infile: {input_path}, outfile: {output_path}")

    # open the json
    with open(input_path, "r") as infile: 
        data = json.load(infile)

    print("opening output csv...")

    # open output csv
    with open(output_path, "w", newline="", encoding="utf-8") as outfile: 

        writer = csv.DictWriter(outfile, fieldnames=list(patterns))

        #write the header
        writer.writeheader() 

        data_list = data.get("data").get("entries")
        print(f"data size: {len(data_list)}")

        for element in data.get("data").get("entries"): 

            body = element.get("body").get("content")

            # extract data for row 
            row = extract_fields(body, patterns) 

            if require_all and any(val is None for val in row.values()): 
                skipped += 1 
                continue 

            writer.writerow(row)
            written += 1 

    return written, skipped

path_infile = "/home/seth-hall/j_research_desktop/APEX_replay/data/apex_end_of_run_logs.json"

if __name__ == "__main__":

    written, skipped = json_to_csv(path_infile, "logs/halog_job_data.csv", require_all=True)

    print(f"Written: {written}, skipped {skipped}")
