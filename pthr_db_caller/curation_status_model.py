from pthr_db_caller.db_caller import DBCaller


class CurationStatus:
    def __init__(self, curation_status_id):
        self.curation_status_id = curation_status_id


class PaintCurationStatusHelper:
    def __init__(self, curation_status_tablename, classification_version_sid):
        self.db_caller = DBCaller()
        self.curation_status_tablename = curation_status_tablename
        self.classification_version_sid = classification_version_sid

    def get_family_classification_id(self, family_id):
        query = """
        select classification_id from panther_upl.classification
        where classification_version_sid = {classification_version_sid}
        and accession = '{family_id}';
        """.format(classification_version_sid=self.classification_version_sid, family_id=family_id)

        results = self.db_caller.run_cmd_line_args(query.rstrip(), no_header_footer=True)
        if len(results[1:]) > 0:
            return results[1][0]
        else:
            raise Exception("No classification found for '{}'".format(family_id))

    def insert_curation_status(self, family_id, status_id):
        family_cls_id = self.get_family_classification_id(family_id)
        # Insert new curation_status

        query = """
        INSERT INTO panther_upl.{curation_status_tablename}
        (curation_status_id, status_type_sid, classification_id, user_id, 
            creation_date)
        VALUES (nextval('uids'), {status_id}, {family_cls_id}, 1113, now())
        """.format(curation_status_tablename=self.curation_status_tablename, status_id=status_id,
                   family_cls_id=family_cls_id)

        self.db_caller.run_cmd_line_args(query.rstrip(), no_header_footer=True)
