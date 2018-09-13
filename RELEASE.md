# Release Instructions

1. Change the `VERSION` file to drop the ".dev0" and get the changed merged into master
2. Create a release in Github using the new version as the title
3. Run the deployment project in Bamboo for the build that includes the version change
4. Bump the version in the `VERSION` file and add the ".dev0" extension
