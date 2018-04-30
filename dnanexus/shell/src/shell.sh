#!/bin/bash

main() {

    echo "Hours to live: '$hours_to_live'"

    sleep $(($hours_to_live * 3600))

}
