import json
import logging
import re
import time
from datetime import datetime

import matplotlib.pyplot as plt
import pandas as pd
import watchdog.events
from matplotlib.backends.backend_pdf import PdfPages
from simple_salesforce import Salesforce
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError
from watchdog.observers.polling import PollingObserver as Observer

### Salesforce API Login ###
logininfo = json.load(open('/ghds/groups/labdesk/bshih/salesforce_login.json'))

username = logininfo['username']
password = logininfo['password']
security_token = logininfo['security_token']
domain = 'login'

sf = Salesforce(username=username, password=password, security_token=security_token)

### Slack Token ###
slack_token = json.load(open('/ghds/groups/labdesk/bshih/slack_login.json'))['SLACK_TOKEN']
client = WebClient(token=slack_token)

### Logger Setup ###
LOGGER_NAME = 'clinical_sample_observation'
LOG_FILE_DIR = '/ghds/groups/labdesk/bshih/clinical_observation'

logger = logging.getLogger(LOGGER_NAME)

log_file = f"{LOG_FILE_DIR}/clinical_sample_observation.log"
handler = logging.FileHandler(log_file)

handler.setLevel(logging.INFO)
logger.setLevel(logging.INFO)

formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", datefmt='%Y-%m-%d %H:%M:%S')
handler.setFormatter(formatter)
logger.addHandler(handler)


class Handler(watchdog.events.FileSystemEventHandler):
    def __init__(self):
        watchdog.events.FileSystemEventHandler()
        logger.info(f'Program has been initiated!')

    @staticmethod
    def sample_finder(all_samples):
        """
        :param all_samples:
        :return: Filtered dataframe with samples whose RNAse <= 50 and G19 >= 0.01, singlicates.
        """
        samples = all_samples[all_samples['run_sample_id'].str.startswith(('G', 'H'), na=False)]
        samples = samples[~samples['run_sample_id'].str.startswith('Ht', na=False)]
        filtered = samples.groupby('run_sample_id').filter(
            lambda x: (x['rnase_count'].median() <= 50) & (x['covid_ratio'].median() >= 0.01))

        return filtered

    @staticmethod
    def salesforce_query(samples):
        """
        :param samples:
        :return: Salesforce query of filtered samples
        """
        weird_samples = []
        for i in samples['run_sample_id'].unique():
            weird_samples.extend(sf.query_all(f"SELECT GH_Sample_ID__c, Status, Specimen_Collection_Date_Time__c, Specimen_receipt_date__c,\
                      State_Authorities_Notified_Date__c, Site_Name__c \
                      FROM Order WHERE GH_Sample_ID__c = '{i}'").get('records'))

        dataframe = pd.DataFrame(weird_samples)
        df = dataframe.drop(columns='attributes').rename(columns={'GH_Sample_ID__c': 'GH Sample ID',
                                                                  'Specimen_Collection_Date_Time__c': 'Specimen Collection Date/Time',
                                                                  'Specimen_receipt_date__c': 'Specimen receipt date',
                                                                  'State_Authorities_Notified_Date__c': 'State Authorities Notified Date',
                                                                  'Site_Name__c': 'Site Name'})
        df['Specimen Collection Date/Time'] = pd.to_datetime(df['Specimen Collection Date/Time']).dt \
            .tz_convert('US/Pacific').dt.strftime('%-m/%-d/%Y, %-I:%M %p')
        df['State Authorities Notified Date'] = pd.to_datetime(df['State Authorities Notified Date']).dt \
            .tz_convert('US/Pacific').dt.strftime('%-m/%-d/%Y, %-I:%M %p')
        df['Specimen receipt date'] = pd.to_datetime(df['Specimen receipt date']).dt.strftime('%-m/%-d/%y')

        return df

    @staticmethod
    def merge_sample_salesforce(samples, salesforce_data):
        """
        :param samples:
        :param salesforce_data:
        :return: Merged dataframe samples with salesforce data
        """
        df = pd.merge(samples, salesforce_data, left_on='run_sample_id', right_on='GH Sample ID', how='left').loc[:,
             ['runid',
              'run_sample_id',
              'covid_ratio',
              'covid_count',
              'rnase_count',
              'spikein_count',
              'replicate_call',
              'replicate_flags',
              'Specimen Collection Date/Time',
              'Specimen receipt date',
              'State Authorities Notified Date',
              'Site Name']] \
            .sort_values(by=['runid', 'run_sample_id', 'covid_ratio'], ascending=[False, True, True])

        median_scores = df.groupby('run_sample_id').median().rename(columns={'covid_ratio': 'new_median'})[
            ['new_median']].round(2)
        df = pd.merge(df, median_scores, left_on='run_sample_id', right_on='run_sample_id', how='left')
        median = df.pop('new_median')
        df.insert(2, 'new_median', median)

        df['group'] = (df['run_sample_id'].shift() != df['run_sample_id']).cumsum()
        df.loc[df.duplicated('group'), ['runid', 'run_sample_id', 'new_median', 'Specimen Collection Date/Time',
                                        'Specimen receipt date', 'State Authorities Notified Date', 'Site Name',
                                        'new_median']] = ''
        df.drop(columns=['group'], inplace=True)

        df = df.astype({'covid_count': 'int', 'rnase_count': 'int', 'spikein_count': 'int'}).round({'covid_ratio': 2})
        return df

    @staticmethod
    def create_pdf(samples, file_name):
        fig, ax = plt.subplots(figsize=(35, len(samples) // 2))
        ax.axis('tight')
        ax.axis('off')
        the_table = ax.table(cellText=samples.values, colLabels=samples.columns, loc='center', rowLoc='right',
                             bbox=[0, 0, 1, 1])

        the_table.auto_set_font_size(False)
        the_table.set_fontsize(12)

        for i in range(0, len(samples.columns)):
            [the_table[(j + k, i)].set_facecolor("#e0e0e0") for j in range(4, len(samples), 6) for k in range(3)]

        the_table.auto_set_column_width(col=list(range(len(samples.columns))))

        pp = PdfPages(f"/ghds/groups/labdesk/bshih/clinical_observation/{file_name}.pdf")
        pp.savefig(fig, bbox_inches='tight', dpi=300)
        pp.close()

    @staticmethod
    def slack_upload(fcid, file_name):
        filepath = f'/ghds/groups/labdesk/bshih/clinical_observation/{file_name}.pdf'

        try:
            result = client.files_upload(
                channels='#g19_sample_observation',
                file=filepath,
                initial_comment=f"{fcid}: These samples have RNAse Count <= 50 and G19 Score >= 0.01")
            assert result["file"]
            logger.info('File Successfully Uploaded to Slack!')
        except SlackApiError as e:
            # You will get a SlackApiError if "ok" is False
            assert e.response["ok"] is False
            assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
            logger.info(f"Got an error: {e.response['error']}")

    def on_created(self, event):
        today = datetime.strftime(datetime.today(), '%y%m%d')
        if bool(re.match(f'^.*?{today}.*?c19_read_counts\.hdr\.tsv$', event.src_path)):

            flowcell = event.src_path[20:51]
            logger.info(f'Flowcell {flowcell} c19_read_counts {event.event_type}!')

            current_read_counts = pd.read_csv(event.src_path, sep='\t')
            low_rnase_high_covid = Handler.sample_finder(current_read_counts)

            if len(low_rnase_high_covid) != 0:
                logger.info(f'There are anomaly samples in: {flowcell}')

                salesforce_data = Handler.salesforce_query(low_rnase_high_covid)
                merged_data = Handler.merge_sample_salesforce(low_rnase_high_covid, salesforce_data)

                filename = flowcell + '_lowRNAse_COVIDpositive'

                Handler.create_pdf(merged_data, filename)
                Handler.slack_upload(flowcell, filename)

                logger.info(f'Salesforce API Usage: {sf.api_usage}')

            else:
                client.chat_postMessage(channel='#g19_sample_observation',
                                        text=f'{flowcell}: No QC fail samples that are covid positive today! Have a great day :)')
                logger.info('No QC fail samples that are covid positive today! Have a great day :)')


if __name__ == "__main__":
    src_path = "/ghds/cv19/analysis"
    event_handler = Handler()
    observer = Observer()
    observer.schedule(event_handler, path=src_path, recursive=True)
    observer.start()
    try:
        while True:
            logger.info("Program is running every 1 Hour!")
            time.sleep(3600)
    except KeyboardInterrupt:
        observer.stop()
        observer.join()
