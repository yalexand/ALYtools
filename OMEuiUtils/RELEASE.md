Release steps:

 * Merge all PRs
 * Create a tag on bintray
 * mvn release:prepare (if you haven't already)
 * mvn release:perform
 * Go to bintray and "publish".

If a release needs to be rolled back:
 * clean up the files created by mvn release:...
 * Revert master
 * Delete the tag locally & remotely
