set -e


print_message() {
  local message="$1"
  local severity="$2"
  local red='\033[0;31m'
  local green='\033[0;32m'
  local yellow='\033[1;33m'
  local nc='\033[0m'

  case "${severity}" in
    "info" ) echo -e "${nc}${message}${nc}";;
    "ok" ) echo -e "${green}${message}${nc}";;
    "error" ) echo -e "${red}${message}${nc}";;
    "warn" ) echo -e "${yellow}${message}${nc}";;
    *) echo -e "${message}";; # Default case if severity is not recognized
  esac
}

# check if release.tar.gz exists
if [ -f release.tar.gz ]; then
    print_message "release.tar.gz already exists" "warn"
    exit 1
fi

# check if wget is installed
if ! [ -x "$(command -v wget)" ]; then
    print_message "wget is not installed" "warn"
    exit 1
fi

# check if wdl-hub directory exists
if [ -d wdl-hub ]; then
    print_message "wdl-hub directory already exists. Remove it to install a new release" "warn"
    exit 1
fi

# witch github release to download? default to v0.1.3 (ex: v0.1.3)
echo "Which release do you want to download? (default: v0.1.3)"
read release
# if empty then default to v0.1.3
if [ -z "$release" ]; then
    release="v0.1.4"
fi

# print the release
echo "Downloading release $release"
wget https://github.com/lmtani/wdl-hub/releases/download/$release/release.tar.gz
tar -xf release.tar.gz
rm release.tar.gz


print_message "Done!" "ok"

print_message "To import tasks from the WDL Hub, add lines like these into you WDL file:" "info"

print_message "import \"wdl-hub/tasks/wget.wdl\"" "ok"

print_message "And then you can call \"wget.Wget\" in your workflow" "info"

print_message "Once you are done, you can build you release with pumbaa build. It will create the .wdl and .zip files" "info"

print_message "Do you want to add the wdl-hub directory to your .gitignore file? (y/n)" "info"
read answer
if [ "$answer" == "y" ]; then
    print_message "Adding wdl-hub to .gitignore" "ok"
    echo "wdl-hub" >> .gitignore
fi

print_message "Suggestion:" "ok"
print_message "- If you need a highly custom task that is not in the WDL Hub, you can create a new task in a \"tasks\" directory and then import it in your WDL file." "info"
print_message "- If you want to contribute to the WDL Hub with a generic task, you can fork the repository and then create a pull request." "info"
