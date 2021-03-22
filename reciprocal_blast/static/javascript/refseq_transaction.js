document.addEventListener('DOMContentLoaded',function () {
    //dataTables: https://datatables.net/
    //DataTables is included via javascript files in the HTML document refseq_database_transactions
    //two features are used to enable an easy interaction with the available refseq files
    //full description of used input is given in the buttons and select Extension section of
    //the DataTables homepage
    var table = $('#myTable').DataTable( {
        dom: 'Bfrtip',
        "lengthMenu": [ 10, 25, 50, 75, 100 ],
        buttons: [
            'copy',
            'csv',
            'selectAll',
            'selectNone',
            'selectRows'
        ],
        select: true
    } );
   document.getElementById('refseqTable').style.display = "block";
})

//Remove all unselected rows from submit
function submitSelectedTableRows(){
    //showLoader('loading_refseq_table','downloadChoices')
    document.getElementById('refseqTable').style.display = "none";
    var loadingDiv = document.getElementById('loading_refseq_table')
    var downloadButton = document.getElementById('downloadChoices')
    downloadButton.style.display = "none";
    loadingDiv.style.display = "block";
    var table = $('#myTable').DataTable();
    var data = table.rows({selected:true}).data().toArray();
    for(var i=0;i < data.length;i++){
        document.getElementById('tableForm').innerHTML += '<input type="text" name="ftp_path[]" value='+data[i][5] +'>'
        //console.log(data[i][5])
    }
    //var tableRows = table.getElementsByTagName('tbody')[0].getElementsByTagName('tr');
    //table.rows({selected:false} ).remove().draw();
    document.getElementById('tableForm').submit()
}

//Not in use:
//This function can be used to loop over table rows
// Add: <input type="text" id="myInput" onkeyup="searchTable()" placeholder="Search for names..">
function searchTable() {
  // Declare variables
  var input, filter, table, tr, td, i, txtValue;
  input = document.getElementById("myInput");
  filter = input.value.toUpperCase();
  table = document.getElementById("myTable");
  tr = table.getElementsByTagName("tr");

  // Loop through all table rows, and hide those who don't match the search query
    for(j = 0; j < table.rows.length;j++){
        //console.log(j)
        for (i = 0; i < tr.length; i++) {

            td = tr[j].getElementsByTagName("td")[i];

            if (td) {
                txtValue = td.textContent || td.innerText;
                if (txtValue.toUpperCase().indexOf(filter) > -1) {
                    tr[j].style.display = "";
                    break;
                } else {
                    tr[j].style.display = "none";
                }
            }
        }

    }
}
//loader screen with giphy
function showLoader(idLoadingDiv,idButton){
    var loadingDiv = document.getElementById(idLoadingDiv)
    var downloadButton = document.getElementById(idButton)
    downloadButton.style.display = "none";
    loadingDiv.style.display = "block";
};