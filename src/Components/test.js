(async () => {
    const fetch = (await import('node-fetch')).default;

    const sequence = "AGAGTCCTGAGCTGAACCAAGAAGGAGGAGGGGGTCGGGCCTCCGAGGAAGGCCTAGCCGCTGCTGCTGC";
    const program = "blastn";
    const database = "nt";

    const submitBlast = async () => {
        const url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi";
        const params = new URLSearchParams({
            CMD: "Put",
            DATABASE: database,
            PROGRAM: program,
            QUERY: sequence,
            FORMAT_TYPE: "JSON2"
        });

        const response = await fetch(url, {
            method: "POST",
            headers: { "Content-Type": "application/x-www-form-urlencoded" },
            body: params.toString()
        });

        const result = await response.text();
        const ridMatch = result.match(/RID = (\w+)/);
        return ridMatch ? ridMatch[1] : null;
    };

    const checkBlastStatus = async (rid) => {
        const url = `https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=${rid}&FORMAT_OBJECT=SearchInfo`;

        let status = "WAITING";
        while (status === "WAITING") {
            const response = await fetch(url);
            const resultText = await response.text();

            if (resultText.includes("Status=WAITING")) {
                console.log("Job is still running...");
                await new Promise(resolve => setTimeout(resolve, 10000));
            } else if (resultText.includes("Status=FAILED")) {
                console.error("BLAST job failed.");
                return null;
            } else if (resultText.includes("Status=READY")) {
                console.log("BLAST job is ready!");
                status = "READY";
            }
        }

        const resultsUrl = `https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=${rid}&FORMAT_TYPE=JSON2`;
        const resultsResponse = await fetch(resultsUrl);

        // Check if the response is JSON
        const contentType = resultsResponse.headers.get("content-type");
        if (contentType && contentType.includes("application/json")) {
            const results = await resultsResponse.json();
            console.log("BLAST Results:", results);
            return results;
        } else {
            const text = await resultsResponse.text();
            console.error("Unexpected response format. Response text:", text);
            return null;
        }
    };

    const rid = await submitBlast();
    if (rid) {
        await checkBlastStatus(rid);
    }
})();
