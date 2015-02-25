<?php

/**
 * @file
 * Map files using the remote_stream_wrapper module.
 */

/**
 * Class for mapping values using the remote_stream_wrapper_module.
 */
class RfsFeedsProvider extends RfsFeedsBaseProvider {
  protected $uri;

  /**
   * The value is the remote URI.
   */
  public function __construct($value, $config = array()) {
    $this->uri = file_stream_wrapper_uri_normalize($value);
  }

  /**
   * Create a file object like remote_stream_wrapper_file_create_by_uri does,
   * only don't save it directly.
   */
  public function getFileObject() {
    $query = parse_url($this->uri, PHP_URL_QUERY);
    if ($query) {
      $query_array = explode('=', $query);
      $name_encoded = end($query_array);
      $filename = urldecode($name_encoded);
    } else {
      $filename = urldecode(basename($this->uri));
    }
    $file = new stdClass();
    $file->fid = NULL;
    $file->uri = $this->uri;
    $file->filename = $filename;
    $file->filemime = file_get_mimetype($file->uri);
    $file->uid = $GLOBALS['user']->uid;
    $file->status = FILE_STATUS_PERMANENT;
    //$file->display = 1;
    return $file;
  }

  /**
   * Save the file or return an existing one.
   */
  public function save() {
    // Make sure no file with that URI exists.
    $existing = remote_stream_wrapper_file_load_by_uri($this->uri);
    if ($existing) {
      return $existing;
    }

    // Save the file.
    $file = $this->getFileObject();
    file_save($file);
    return $file;
  }
  
}
