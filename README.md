# eXFrame
eXFrame is a reusable framework for building genomics experiments repositories.

## CONTENTS OF THIS FILE
---------------------

  * [Introduction](introduction)
  * [Requirements](requirements)
  * [Recommended Modules](recommended-modules)
  * [Installation](installation)
  * [Configuration](configuration)
  * [Support](support)
  * [Contributing](contributing)
  * [Maintainers](maintainers)
  * [License](license)
  * [Citations](citations)
  
### INTRODUCTION
------------

  * For a full description, please visit:
    * http://www.jbiomedsem.com/content/5/S1/S3
    * http://www.biomedcentral.com/1471-2105/12/452
  * To submit bug reports or feature suggestions:
    * http://github.com/mindinformatics/exframe/issues

### REQUIREMENTS
------------

  * R (&gt;= {Version?}) - https://www.r-project.org/
  * Web Server (Apache,Nginx) 
  * Database Server (MySQL / MariaDB {Version?}, postgreSQL {???}, SQLite {???})
  * PHP (&gt;= 5.2)
  * Drush - http://www.drush.org

### RECOMMENDED MODULES
-------------------


### INSTALLATION
------------

  1. Verify installation dependencies have been installed (R, Drush, Database Server, PHP)
  2. https://www.drupal.org/documentation/install
  3. Log in
  4. Enable eXFrame modules (Administer -> Modules)
    * Check the box next to each module that starts with "exframe_". Click "Save Configuration" at the bottom of the page.
  5. Create Content
    * Import Taxonomies
    * Create Users
      * Administer -> People -> Add Users
      * (Or import CSV)
    * Create Group(s)
      * Administer -> Content -> Add Content -> Group
    * Create Project(s)
      * Administer -> Content -> Add Content -> Project
    * Create Experiment(s)
      * Administer -> Content -> Add Content -> Experiment
    * Create Sample(s)
      * Administer -> Content -> Add Content -> Sample
    * Create Sample Group(s)
      * Administer -> Content -> Add Content -> Sample Group
  
### CONFIGURATION
-------------


### SUPPORT
-------


#### TROUBLESHOOTING
---------------


#### FAQ
---

  * Q: 
    * A: 


### CONTRIBUTING
------------

#### ROADMAP
------------

  * [] Generate "install profile"
  * [] Update code (style) to meet Drupal coding standards
  * [] Create installation script (drush)
  * [] Create demo installation (eXFrame example module / sub-profile)
  * [] Improve documentation
    * [] Installation
    * [] Configuration
    * [] Customization
  * [] Continuous Integration (Jenkins, Travis)
    * [] Regression Testing
  * [] Create "Official" Drupal project
  
### LICENSE
-------
This is licensed under GPL v2
See [LICENSE.txt](./LICENSE.txt) for the full text

### MAINTAINERS (CONTRIBUTORS)
-----------

Current Maintainers:
  * First M. Last (uname) - webpage_link
  * ...
This project has been sponsored by / funded through:
  * ...

### CITATION
--------

Sinha, A., Merrill, E., Armstrong, S., Clark, T., & Das, S. (2011). EXframe: Reusable framework for storage, analysis and visualization of genomics experiments. BMC Bioinformatics, 12, 452-452. doi:10.1186/1471-2105-12-452

http://www.biomedcentral.com/1471-2105/12/452


Merrill, E., Corlosquet, S., Ciccarese, P., Clark, T., & Das, S. (2014). Semantic Web repositories for genomics data using the eXframe platform. Journal of Biomedical Semantics, 5(Supplement 1), S3-S3. doi:10.1186/2041-1480-5-S1-S3 

http://www.jbiomedsem.com/content/5/S1/S3
