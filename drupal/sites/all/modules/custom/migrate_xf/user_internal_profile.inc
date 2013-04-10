<?php
/**
 * Migration class to import user profiles.
 */
class ExframeInternalProfileMigration extends Migration {
  public function __construct() {
    parent::__construct();
    $this->description = t('Internal user profiles import from CSV');
    $this->dependencies = array('ExframeUser');

    $this->map = new MigrateSQLMap($this->machineName,
      array(
        'uid' => array(
          'type' => 'int',
          'unsigned' => TRUE,
          'not null' => TRUE,
        )
      ),
      MigrateDestinationProfile2::getKeySchema()
    );

    // Our data is in a CSV file.
    $this->source = new MigrateSourceCSV(DRUPAL_ROOT . '/../migration/data/profiles_internal.csv', array(), array('header_rows' => 1), $this->fields());
    $this->destination = new MigrateDestinationProfile2('xf_profile');

    // Core mappings.

    // Field mappings.
    $this->addFieldMapping('uid', 'uid');
    $this->addFieldMapping('field_xfp_firstname', 'FirstName');
    $this->addFieldMapping('field_xfp_minitial', 'MInitial');
    $this->addFieldMapping('field_xfp_lastname', 'LastName');
    $this->addFieldMapping('field_xfp_affiliation', 'Org');
    $this->addFieldMapping('field_xfp_department', 'Department');
    $this->addFieldMapping('field_xfp_lab', 'Lab');
    $this->addFieldMapping('field_xfp_phone', 'Phone');

    $arguments = array(
      'thoroughfare' => array('source_field' => 'Street1'),
      'premise' => array('source_field' => 'Street2'),
      'locality' => array('source_field' => 'City'),
      'postal_code' => array('source_field' => 'Postal_Code'),
      'administrative_area' => array('source_field' => 'State'),
    );
    // The country should be passed in as the primary value.
    $this->addFieldMapping('field_xfp_address', 'profile_country')
         ->arguments($arguments)
         ->defaultValue('US');
  }

  // Optional - keys are field names, values are descriptions. Use to override
  // the default descriptions, or to add additional source fields which the
  // migration will add via other means (e.g., prepareRow()).
  function fields() {
    return array(
      'uid' => 'User ID',
    );
  }
}
