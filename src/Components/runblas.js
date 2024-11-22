const axios = require('axios');

const sequence = "AGAGTCCTGAGCTGAACCAAGAAGGAGGAGGGGGTCGGGCCTCCGAGGAAGGCCTAGCCGCTGCTGCTGC"; // Replace with the actual DNA sequence

async function queryBlast() {
    try {
        // Step 1: Submit the sequence for BLAST search
        let submitResponse = await axios.post('https://blast.ncbi.nlm.nih.gov/Blast.cgi', null, {
            params: {
                CMD: 'Put',
                DATABASE: 'nt',
                PROGRAM: 'blastn',
                QUERY: sequence,
                FORMAT_TYPE: 'XML'
            }
        });

        // Extract the Request ID (RID) and Estimated Time (RTOE)
        const rid = submitResponse.data.match(/RID = (\w+)/)[1];
        const rtoe = parseInt(submitResponse.data.match(/RTOE = (\d+)/)[1]);

        console.log(`Request ID (RID): ${rid}`);
        console.log(`Estimated Time to Completion (RTOE): ${rtoe} seconds`);

        // Step 2: Check the status of the request
        let status = 'WAITING';
        while (status === 'WAITING') {
            await new Promise(resolve => setTimeout(resolve, rtoe * 1000)); // Wait based on RTOE
            let statusResponse = await axios.get('https://blast.ncbi.nlm.nih.gov/Blast.cgi', {
                params: {
                    CMD: 'Get',
                    FORMAT_OBJECT: 'SearchInfo',
                    RID: rid
                }
            });

            if (statusResponse.data.includes("Status=WAITING")) {
                console.log("Job is still processing...");
            } else if (statusResponse.data.includes("Status=FAILED")) {
                console.error("BLAST search failed.");
                return;
            } else {
                status = 'READY';
            }
        }

        // Step 3: Retrieve the results
        let resultResponse = await axios.get('https://blast.ncbi.nlm.nih.gov/Blast.cgi', {
            params: {
                CMD: 'Get',
                FORMAT_TYPE: 'XML',
                RID: rid
            }
        });

        console.log("BLAST results:", resultResponse.data);

    } catch (error) {
        console.error("Error querying BLAST:", error);
    }
}

queryBlast();
