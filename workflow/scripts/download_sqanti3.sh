#!env bash
set -o pipefail

sq3_version="5.3.6"
mkdir -p resources

pushd resources
wget -q -O sq3_tmp.zip https://github.com/ConesaLab/SQANTI3/archive/refs/tags/v${sq3_version}.zip
unzip -q -o sq3_tmp.zip
rm sq3_tmp.zip
rm -rf Sqanti3
mv SQANTI3-${sq3_version} Sqanti3
rm -rf Sqanti3/test
rm -rf Sqanti3/example
rm -rf Sqanti3/data
popd


# wget -q -O sq3_tmp.zip https://github.com/ConesaLab/SQANTI3/releases/download/v5.3.6/SQANTI3_v5.3.6.zip &&
#     unzip -q -o sq3_tmp.zip -d Sqanti3_v5.3.6 &&
#     rm sq3_tmp.zip &&
#     rm -rf Sqanti3_v5.3.6/examplecd