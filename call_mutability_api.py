#import requests

#MUTAGENE_URL = "http://localhost:5000"
#MUTAGENE_URL = "https://www.ncbi.nlm.nih.gov/research/mutagene"

import requests
import sys
MUTAGENE_URL = "https://dev.ncbi.nlm.nih.gov/research/mutagene"
#MUTAGENE_URL = "https://www.ncbi.nlm.nih.gov/research/mutagene"


def retrieve_mutability(gene, model,download):
    """
    """

    url = MUTAGENE_URL + '/api/gene'
    r = requests.post(url, data={
        'gene': gene,
        'signature': str(model),
        'download': str(download),
        'mutations': 'profile',
        'mutation_rate': '0',
        # 'mutation_rate_value': 120,
    })
    if r.status_code == 200:
        if download == "4":
            muts = r.json()
            for entry in muts['data']:
                print("{}\t{}".format(entry['mutation'], entry['impact']))

        else:
            return(r.text)


if __name__ == '__main__':
    gene = sys.argv[1]
    model = sys.argv[2]
    download = sys.argv[3]
    print(retrieve_mutability(gene, model,download))
