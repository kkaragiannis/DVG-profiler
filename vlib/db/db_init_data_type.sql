/*
 *  ::718604!
 * 
 * Copyright(C) November 20, 2014 U.S. Food and Drug Administration
 * Authors: Dr. Vahan Simonyan (1), Dr. Raja Mazumder (2), et al
 * Affiliation: Food and Drug Administration (1), George Washington University (2)
 * 
 * All rights Reserved.
 * 
 * The MIT License (MIT)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

START TRANSACTION;

source db_init_data_include.sql;

SET @type_domain = 1954115685;
SET @type_type_id = 1;

DELETE FROM UPPerm WHERE objID = @type_type_id AND domainID = @type_domain;
DELETE FROM UPObjField WHERE objID = @type_type_id AND domainID = @type_domain;
DELETE FROM UPObj WHERE objID = @type_type_id AND domainID = @type_domain;

CALL sp_obj_create_v2(@system_group_id, @system_membership, @type_domain, @type_type_id, @type_domain, @type_type_id, @system_permission, @system_flags);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'name\',NULL,\'type\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'title\',NULL,\'object type\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'description\',NULL,\'Object type\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'prefetch\',NULL,\'1\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.1.1\',\'created\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.1.1\',\'Created\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.1.1\',\'datetime\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_readonly_fg\',\'1.1.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.1.1\',\'1\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.2.1\',\'description\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.2.1\',\'Description\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.2.1\',\'text\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.2.1\',\'3\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.2.1\',\'1\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.3.1\',\'fields\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.3.1\',\'Fields\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.3.1\',\'array\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.3.1\',\'8\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.3.1\',\'1\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.4.1\',\'field_basics\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.4.1\',\'Basics\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.4.1\',\'list\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.4.1\',\'fields\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.4.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_multi_fg\',\'1.4.1\',\'1\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.5.1\',\'field_brief\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.5.1\',\'Brief format\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.5.1\',\'string\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.5.1\',\'field_presentation\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.5.1\',\'31\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.5.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_description\',\'1.5.1\',\'Empty or 0 not a brief member; [text]$_(v)[text], where $_(v) is substituted by prop value\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.6.1\',\'field_constraint\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.6.1\',\'Constraint type\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.6.1\',\'string\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.6.1\',\'field_constraint_list\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.6.1\',\'21\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.6.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_constraint\',\'1.6.1\',\'choice\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_constraint_data\',\'1.6.1\',\'choice|choice+|regexp|range|url|search|search+|one-of|type|eval\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.7.1\',\'field_constraint_data\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.7.1\',\'Constraint data\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.7.1\',\'text\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.7.1\',\'field_constraint_list\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.7.1\',\'22\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.7.1\',\'1\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.8.1\',\'field_constraint_description\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.8.1\',\'Constraint description\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.8.1\',\'text\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.8.1\',\'field_constraint_list\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.8.1\',\'23\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.8.1\',\'1\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.9.1\',\'field_constraint_list\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.9.1\',\'Value constraints\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.9.1\',\'list\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.9.1\',\'fields\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.9.1\',\'3\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_multi_fg\',\'1.9.1\',\'1\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.10.1\',\'field_default_encoding\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.10.1\',\'Default encoding\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.10.1\',\'integer\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.10.1\',\'field_basics\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.10.1\',\'7\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.10.1\',\'1\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.11.1\',\'field_default_value\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.11.1\',\'Default value\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.11.1\',\'string\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.11.1\',\'field_basics\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.11.1\',\'6\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.11.1\',\'1\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.12.1\',\'field_description\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.12.1\',\'Description\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.12.1\',\'text\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.12.1\',\'field_presentation\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.12.1\',\'33\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.12.1\',\'1\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.13.1\',\'field_flags\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.13.1\',\'Flags\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.13.1\',\'list\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.13.1\',\'fields\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.13.1\',\'2\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_multi_fg\',\'1.13.1\',\'1\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.14.1\',\'field_include_type\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.14.1\',\'Inclusion from type\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.14.1\',\'obj\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.14.1\',\'field_basics\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.14.1\',\'8\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.14.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_constraint\',\'1.14.1\',\'search\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_constraint_data\',\'1.14.1\',\'{\\n  \\\"explorer\\\": true,\\n  \\\"fetch\\\": \\\"id\\\",\\n  \\\"inline\\\": \\\"name\\\",\\n  \\\"inline_url\\\": \\\"http://?cmd=propget&type=^type$&mode=csv&raw=1&useTypeDomainId=1&useTypeUPObj=1\\\",\\n  \\\"url\\\": \\\"http://?cmd=objList&type=^type$&mode=csv&raw=1&useTypeDomainId=1&useTypeUPObj=1\\\"\\n}\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_description\',\'1.14.1\',\'type whose fields will be included into this one as sub-fields of array or list\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.15.1\',\'field_is_batch_fg\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.15.1\',\'Batchable\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.15.1\',\'bool\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.15.1\',\'field_flags\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.15.1\',\'18\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_default_value\',\'1.15.1\',\'0\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.15.1\',\'1\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.16.1\',\'field_is_hidden_fg\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.16.1\',\'Hidden\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.16.1\',\'bool\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.16.1\',\'field_flags\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.16.1\',\'16\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_default_value\',\'1.16.1\',\'0\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.16.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_description\',\'1.16.1\',\'field is not visible on pages\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.17.1\',\'field_is_key_fg\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.17.1\',\'Key\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.17.1\',\'bool\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.17.1\',\'field_flags\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.17.1\',\'12\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_default_value\',\'1.17.1\',\'0\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.17.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_description\',\'1.17.1\',\'field is part of unique identifier of the record or parent list or array\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.18.1\',\'field_is_multi_fg\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.18.1\',\'Multi-value\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.18.1\',\'bool\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.18.1\',\'field_flags\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.18.1\',\'15\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_default_value\',\'1.18.1\',\'0\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.18.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_description\',\'1.18.1\',\'field can have multiple values\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.19.1\',\'field_is_optional_fg\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.19.1\',\'Optional\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.19.1\',\'bool\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.19.1\',\'field_flags\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.19.1\',\'14\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_default_value\',\'1.19.1\',\'0\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.19.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_description\',\'1.19.1\',\'field could be missing completely\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.20.1\',\'field_is_readonly_fg\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.20.1\',\'Read-only\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.20.1\',\'integer\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.20.1\',\'field_flags\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.20.1\',\'13\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_default_value\',\'1.20.1\',\'0\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.20.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_constraint\',\'1.20.1\',\'choice\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_constraint_data\',\'1.20.1\',\'0///normal field|1///calculated outside and NOT editable on web pages|-1///write-once|-2///submit-once|2///autofill with default value\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.21.1\',\'field_is_summary_fg\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.21.1\',\'Summary\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.21.1\',\'bool\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.21.1\',\'field_flags\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.21.1\',\'17\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_default_value\',\'1.21.1\',\'0\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.21.1\',\'1\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.22.1\',\'field_is_virtual_fg\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.22.1\',\'Virtual\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.22.1\',\'bool\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.22.1\',\'field_flags\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.22.1\',\'15\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_default_value\',\'1.22.1\',\'0\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.22.1\',\'1\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.23.1\',\'field_link_url\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.23.1\',\'Link URL\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.23.1\',\'string\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.23.1\',\'field_presentation\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.23.1\',\'32\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.23.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_description\',\'1.23.1\',\'next to the field, adds a link using current field value as $(v) substitution in url\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.24.1\',\'field_name\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.24.1\',\'Name\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.24.1\',\'string\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.24.1\',\'field_basics\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.24.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_key_fg\',\'1.24.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_description\',\'1.24.1\',\'WARNING: changing the name in an existing type will break existing objects of this type or descendent types!!! Case-insensitive unique name within this type, [a-zA-Z0-9_-] character set, must not be used by any other field within this type, must not start with _\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.25.1\',\'field_order\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.25.1\',\'Order\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.25.1\',\'real\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.25.1\',\'field_basics\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.25.1\',\'5\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.25.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_description\',\'1.25.1\',\'numeric (integer or real) order for sorting within parent field\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.26.1\',\'field_parent\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.26.1\',\'Parent field\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.26.1\',\'string\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.26.1\',\'field_basics\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.26.1\',\'4\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.26.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_description\',\'1.26.1\',\'name of parent field\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.27.1\',\'field_presentation\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.27.1\',\'Presentation\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.27.1\',\'list\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.27.1\',\'fields\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.27.1\',\'4\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_multi_fg\',\'1.27.1\',\'1\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.28.1\',\'field_role\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.28.1\',\'Role\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.28.1\',\'string\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.28.1\',\'field_flags\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.28.1\',\'11\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.28.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_constraint\',\'1.28.1\',\'choice\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_constraint_data\',\'1.28.1\',\'in|out\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.29.1\',\'field_title\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.29.1\',\'Title\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.29.1\',\'string\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.29.1\',\'field_basics\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.29.1\',\'2\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.29.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_description\',\'1.29.1\',\'human-readable field title to show on pages\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.30.1\',\'field_type\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.30.1\',\'Value type\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.30.1\',\'string\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_parent\',\'1.30.1\',\'field_basics\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.30.1\',\'3\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_constraint\',\'1.30.1\',\'choice\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_constraint_data\',\'1.30.1\',\'string|text|integer|real|bool|array|list|date|time|datetime|url|obj|password|file|type2array///Include an object type as array|type2list///Include an object type as list\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.31.1\',\'is_abstract_fg\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.31.1\',\'Abstract\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.31.1\',\'bool\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.31.1\',\'5\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_default_value\',\'1.31.1\',\'0\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.31.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_description\',\'1.31.1\',\'if a type is abstract, objects of that type cannot exist in the system - but such a type can be used for multiple inheritance\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.32.1\',\'modified\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.32.1\',\'Modified\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.32.1\',\'datetime\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_readonly_fg\',\'1.32.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.32.1\',\'1\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.33.1\',\'name\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.33.1\',\'Name\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.33.1\',\'string\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.33.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_key_fg\',\'1.33.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_summary_fg\',\'1.33.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_brief\',\'1.33.1\',\'$_(v) <b>Type</b>\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_description\',\'1.33.1\',\'case-insensitive unique name, [a-zA-Z0-9_-] character set, must not be used by any other type, must not start with _\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.34.1\',\'parent\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.34.1\',\'Parent type(s)\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.34.1\',\'obj\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.34.1\',\'4\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.34.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_multi_fg\',\'1.34.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_constraint\',\'1.34.1\',\'search\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_constraint_data\',\'1.34.1\',\'{\\n  \\\"explorer\\\": true,\\n  \\\"fetch\\\": \\\"id\\\",\\n  \\\"inline\\\": \\\"name\\\",\\n  \\\"inline_url\\\": \\\"http://?cmd=propget&type=^type$&mode=csv&raw=1&useTypeDomainId=1&useTypeUPObj=1\\\",\\n  \\\"url\\\": \\\"http://?cmd=objList&type=^type$&mode=csv&raw=1&useTypeDomainId=1&useTypeUPObj=1\\\"\\n}\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_description\',\'1.34.1\',\'one or more types whose fields will be inherited by this one; only one of them may be non-abstract\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.35.1\',\'prefetch\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.35.1\',\'Pre-loaded\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.35.1\',\'bool\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.35.1\',\'7\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_default_value\',\'1.35.1\',\'0\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_is_optional_fg\',\'1.35.1\',\'1\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_description\',\'1.35.1\',\'should this type and its fields always be loaded by all backend services?\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.36.1\',\'singleton\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.36.1\',\'Singleton\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.36.1\',\'integer\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.36.1\',\'6\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_default_value\',\'1.36.1\',\'0\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_constraint\',\'1.36.1\',\'choice\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_constraint_data\',\'1.36.1\',\'0///No|1///One per owner|2///One system wide\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_description\',\'1.36.1\',\'is only one object of this type allowed to exist?\',NULL,NULL)'), NULL);

CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_name\',\'1.37.1\',\'title\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_title\',\'1.37.1\',\'Title\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_type\',\'1.37.1\',\'string\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_order\',\'1.37.1\',\'2\',NULL,NULL)'), NULL);
CALL sp_obj_prop_set_v3(@system_group_id, @system_membership, @type_domain, @type_type_id, @system_permission, CONCAT('(', @type_domain, ',', @type_type_id, ',\'field_description\',\'1.37.1\',\'human-readable title to show on pages\',NULL,NULL)'), NULL);


COMMIT;
