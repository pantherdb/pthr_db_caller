from pthr_db_caller.db_caller import DBCaller


class Comment:
    def __init__(self, comment_id, classification_id, remark, created_by, creation_date, obsoleted_by,
                 obsolescence_date, protein_id=None, node_id=None):
        self.comment_id = comment_id
        self.classification_id = classification_id
        self.protein_id = protein_id
        self.remark = remark
        self.created_by = created_by
        self.creation_date = creation_date
        self.obsoleted_by = obsoleted_by
        self.obsolescence_date = obsolescence_date
        self.node_id = node_id


def parse_results_to_comments(results):
    comments = []
    for r in results:
        comment = Comment(
            comment_id=r[0],
            classification_id=r[1],
            protein_id=r[2],
            remark=r[3],
            created_by=r[4],
            creation_date=r[5],
            obsoleted_by=r[6],
            obsolescence_date=r[7],
            node_id=r[8]
        )
        comments.append(comment)
    return comments


def parse_results_to_fam_data_struct(results):
    # query field list must be: family, n.public_id, gc.accession, cc.confidence_code, q.qualifier

    # Data model
    # families = {
    #     "PTHR10000": {
    #         "PTN000389662": [
    #             {
    #                 "term": "GO:0060090",
    #                 "ev_code": "IBD"
    #             },
    #             {
    #                 "term": "GO:0010739",
    #                 "ev_code": "IBD"
    #             }
    #         ]
    #     }
    # }
    families = {}
    for r in results:
        family = r[0]
        node_ptn = r[1]
        term = r[2]
        ev_code = r[3]
        qualifier = r[4]
        if qualifier is None:
            qualifier = ""
        annotation = {"term": term, "ev_code": ev_code, "qualifier": qualifier}
        if family not in families:
            families[family] = {}
        if node_ptn not in families[family]:
            families[family][node_ptn] = [annotation]
        else:
            families[family][node_ptn].append(annotation)
    return families


class PthrCommentHelper:
    def __init__(self, comments_tablename, classification_version_sid):
        self.db_caller = DBCaller()
        self.comments_tablename = comments_tablename
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

    def get_comments(self, family_id):
        # Return record
        query = """
        select cm.* from panther_upl.{comments_tablename} cm
        join panther_upl.classification c on c.classification_id = cm.classification_id
        where c.classification_version_sid = {classification_version_sid}
        and c.accession = '{family_id}';
        """.format(comments_tablename=self.comments_tablename,
                   classification_version_sid=self.classification_version_sid, family_id=family_id)

        # print(query)
        results = self.db_caller.run_cmd_line_args(query.rstrip(), no_header_footer=True)
        return parse_results_to_comments(results[1:])

    def update_comment(self, family_cls_id, comment_text):
        # Append to remark field
        query = """
        update panther_upl.{comments_tablename} cm
        set remark = remark || '\n' || current_date || ': {comment_text}\n'
        where cm.classification_id = {family_cls_id};
        """.format(comments_tablename=self.comments_tablename, comment_text=comment_text, family_cls_id=family_cls_id)

        self.db_caller.run_cmd_line_args(query.rstrip(), no_header_footer=True)

    def insert_comment(self, family_cls_id, comment_text):
        # Insert new comment
        query = """
        insert into panther_upl.{comments_tablename}
        (comment_id, classification_id, protein_id, remark, created_by, 
            creation_date, obsoleted_by, obsolescence_date, node_id)
        VALUES (nextval('uids'), {family_cls_id}, null, current_date || ': {comment_text}\n', 1113, now(), null, null, null);
        """.format(comments_tablename=self.comments_tablename, family_cls_id=family_cls_id, comment_text=comment_text)

        self.db_caller.run_cmd_line_args(query.rstrip(), no_header_footer=True)

    def update_or_insert_comment(self, family_id, comment_text):
        existing_comments = self.get_comments(family_id)
        if len(existing_comments) > 0:
            # Can update
            cls_id = existing_comments[0].classification_id
            self.update_comment(cls_id, comment_text)
        else:
            # Gonna have to get family classification_id
            cls_id = self.get_family_classification_id(family_id)
            self.insert_comment(cls_id, comment_text)

    # TODO: Standardize annotation comment lines from annotation query
    # def generate_comments_by_query(self, query, comment_header):
    #     families = parse_results_to_fam_data_struct(results[1:])

# Testing
# helper = PthrCommentHelper("comments", 26)
# test_comments = helper.get_comments("PTHR10202")
