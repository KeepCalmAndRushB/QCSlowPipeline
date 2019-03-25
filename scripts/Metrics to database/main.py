import pandas as pd
import mysql.connector
from sqlalchemy import create_engine

def main(qcresults):
    qc = pd.read_csv(qcresults, sep='\t')
    engine = create_engine('mysql+mysqlconnector://root:qx1qx2@127.0.0.1:3306/msqc_slow', echo=False)
    qc.to_sql(name='qc_all_together', con=engine, if_exists='append', index=False)

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description="Extract QC Metrics.")
    parser.add_argument('--qcresults', default='O:/20190319_QX4_MePh_MA_HeLa_500ng_LC11/QC_Results.tab', help="The table with QC results")

    args = parser.parse_args()

    main(qcresults=args.qcresults)