import psi4
import json

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Run psi4 to calculate bond orders')
    parser.add_argument('-i', '--infile', type=str,
                        help='Input JSON file')

    args = parser.parse_args()
    infile = args.infile
    json_data = json.load(open(infile, 'r'))
    
    psi4.json_wrapper.run_json(json_data)
    
    outfile_name = infile.replace('input', 'output')
    outfile = open(outfile_name, 'w')
    json.dump(json_data, outfile, indent=4, sort_keys=True)

    # Write out raw output
    text_output_name = outfile_name.replace('json', 'txt')
    text_outfile = open(text_output_name, 'w')
    text_outfile.write(json_data['raw_output'])
    text_outfile.close()
